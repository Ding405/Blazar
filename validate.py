"""Synthetic environment validation; not a real-data HST fitting pipeline."""
from pathlib import Path
import json
import tempfile
import numpy as np
from astropy.io import fits
from astropy.modeling.models import Sersic2D
from photutils.aperture import CircularAperture, aperture_photometry
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares
import matplotlib


def in_notebook():
    try:
        return get_ipython().kernel is not None
    except (NameError, AttributeError):
        return False


# Preserve the notebook backend; use a headless backend for standalone runs.
if not in_notebook():
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from PIL import Image

y, x = np.mgrid[:65, :65]
psf = np.zeros((65, 65)); psf[32, 32] = 1
psf = gaussian_filter(psf, 1.6)


def model(p):
    host = Sersic2D(amplitude=1, r_eff=p[2], n=2, x_0=32,
                    y_0=32, ellip=0.2, theta=0.3)(x, y)
    host = gaussian_filter(host, 1.6)
    # p[1] is flux within this finite image, not infinite Sersic total flux.
    host /= host.sum()
    return p[0] * psf + p[1] * host + p[3]


def check_fit(parameters, residual):
    if not np.isfinite(parameters).all() or not np.isfinite(residual).all():
        raise ValueError('Nonfinite fit parameters or residuals')
    if not np.allclose(parameters, truth, rtol=0.05):
        raise ValueError(f'Parameter recovery exceeds 5%, including sky: {parameters}')
    reduced_chi2 = float(np.sum((residual / noise_sigma)**2) /
                         (residual.size - len(parameters)))
    if abs(float(residual.mean())) >= 0.05 or not 0.8 < reduced_chi2 < 1.2:
        raise ValueError('Residual statistics inconsistent with synthetic noise')
    return reduced_chi2


truth = np.array([20000., 30000., 8., 10.])
noise_sigma = 0.5
rng = np.random.default_rng(42)
data = model(truth) + rng.normal(0, noise_sigma, x.shape)
fit = least_squares(lambda p: ((model(p)-data)/noise_sigma).ravel(),
                    [18000, 25000, 6, 9],
                    bounds=([0, 0, 1, 0], [1e6, 1e6, 25, 100]))
if not fit.success:
    raise RuntimeError(fit.message)
best_model = model(fit.x)
residual = data - best_model
reduced_chi2 = check_fit(fit.x, residual)

# In a pasted notebook cell __file__ is absent; write relative to its cwd.
base = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
output_root = base / 'outputs'
output_root.mkdir(parents=True, exist_ok=True)
output_dir = Path(tempfile.mkdtemp(prefix='validation-', dir=output_root))
path = output_dir / 'synthetic.fits'
fits.writeto(path, data, checksum=True)
with fits.open(path, checksum=True) as hdul:
    np.testing.assert_array_equal(hdul[0].data, data)
    if hdul[0].verify_checksum() != 1 or hdul[0].verify_datasum() != 1:
        raise ValueError('FITS checksum validation failed')
table = aperture_photometry(data-fit.x[3], CircularAperture((32,32), r=15))
if float(table['aperture_sum'][0]) <= 20000:
    raise ValueError('Unexpected synthetic aperture photometry')

# Save both the individual residual and a comparison figure.
fig, ax = plt.subplots()
v = float(np.max(np.abs(residual)))
im = ax.imshow(residual, origin='lower', cmap='RdBu_r', vmin=-v, vmax=v)
ax.set_title('Synthetic residual'); fig.colorbar(im, ax=ax)
fig.savefig(output_dir / 'residual.png', dpi=120)
plt.close(fig)
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), layout='constrained')
norm = AsinhNorm(linear_width=3, vmin=8, vmax=max(data.max(), best_model.max()))
for ax, arr, title in zip(axes[:2], [data, best_model],
                          ['Simulated observation', 'Best-fit model']):
    ax.imshow(arr, origin='lower', norm=norm, cmap='magma'); ax.set_title(title)
im = axes[2].imshow(residual, origin='lower', cmap='RdBu_r', vmin=-v, vmax=v)
axes[2].set_title('Residual'); fig.colorbar(im, ax=axes[2], shrink=0.8)
for ax in axes:
    ax.set_xlabel('x [pixel]'); ax.set_ylabel('y [pixel]')
fig.suptitle('Synthetic validation only - not real HST data')
comparison = output_dir / 'comparison.png'
fig.savefig(comparison, dpi=120)
plt.close(fig)
for name in ['residual.png', 'comparison.png']:
    with Image.open(output_dir / name) as image:
        image.verify()
    with Image.open(output_dir / name) as image:
        if np.asarray(image.convert('RGB')).std() == 0:
            raise ValueError(f'Blank output: {name}')
summary = {'parameters': fit.x.tolist(),
           'parameter_names': ['nuclear_image_flux', 'host_image_flux', 'r_eff_pixels', 'sky'],
           'truth': truth.tolist(), 'reduced_chi2': reduced_chi2,
           'residual_mean': float(residual.mean()),
           'residual_std': float(residual.std()),
           'scope': 'Synthetic smoke check; fixed PSF and galaxy shape; no scientific uncertainty estimate'}
(output_dir / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print('PASS: parameter recovery including sky, residuals, FITS checksums, photometry, PNG readback')
print('Recovered [nuclear image flux, host image flux, r_eff, sky]:', fit.x)
print(f'Reduced chi-square: {reduced_chi2:.4f}')
print('Outputs retained:', output_dir)
if in_notebook():
    from IPython.display import Image as NotebookImage, display
    display(NotebookImage(filename=str(comparison)))
print('Synthetic check only; not validation of real HST data or parameter uncertainties.')
