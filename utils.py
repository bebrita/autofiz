# Стандартные библиотеки
from typing import Tuple

# Сторонние библиотеки
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches

# Astropy
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.visualization import simple_norm
from scipy.interpolate import interp1d

# Scipy
from scipy.ndimage import (
    gaussian_filter,
    generate_binary_structure,
    binary_erosion,
    label,
    find_objects
)

from scipy.optimize import leastsq, brentq


def visualize_fits(file_path: str) -> None:
    with fits.open(file_path) as hdul:
        print("\nСтруктура FITS файла:")
        data = hdul[0].data
        header = hdul[0].header

        print("\nВажные параметры из заголовка:")
        print(f"Размер изображения: {data.shape}")
        print(f"Координаты (RA, DEC): {header.get('RA', 'не указано')}, {header.get('DEC', 'не указано')}")
        print(f"Дата наблюдения: {header.get('DATE-OBS', 'не указано')}")

        plt.figure(figsize=(12, 6))
        norm = simple_norm(data, 'sqrt', percent=99)
        plt.subplot(1, 2, 1)
        plt.imshow(data, cmap='gray', origin='lower', norm=norm)
        plt.colorbar(label='Интенсивность')
        plt.title(file_path)
        plt.tight_layout()
        plt.show()


def intetn_fits(file_path: str) -> None:
    with fits.open(file_path) as hdul:
        data = hdul[0].data

        plt.figure(figsize=(12, 6))
        plt.subplot(1, 2, 2)
        plt.hist(data.flatten(), bins=100, log=True)
        plt.xlabel('Значение интенсивности')
        plt.ylabel('Частота (логарифмическая)')
        plt.title('Распределение интенсивностей')

        plt.tight_layout()
        plt.show()

        print("\nАнализ интенсивностей:")
        print(f"Минимальная интенсивность: {np.min(data)}")
        print(f"Максимальная интенсивность: {np.max(data)}")
        print(f"Средняя интенсивность: {np.mean(data)}")
        print(f"Медианная интенсивность: {np.median(data)}")
        print(f"Стандартное отклонение интенсивностей: {np.std(data)}")



def crop_fits(input_path: str, output_path: str,
              x_slice: Tuple[int, int],
              y_slice: Tuple[int, int],
              update_wcs: bool = True) -> None:

    with fits.open(input_path) as hdul:
        data = hdul[0].data
        header = hdul[0].header.copy()
        if x_slice[1] <= x_slice[0] or y_slice[1] <= y_slice[0]:
            raise ValueError("Некорректные границы обрезки (end должно быть больше start)")

        if x_slice[1] > data.shape[1] or y_slice[1] > data.shape[0]:
            raise ValueError("Границы обрезки превышают размеры изображения")

        cropped_data = data[y_slice[0]:y_slice[1], x_slice[0]:x_slice[1]]

        if update_wcs and 'CRPIX1' in header and 'CRPIX2' in header:
            header['CRPIX1'] -= x_slice[0]
            header['CRPIX2'] -= y_slice[0]
            header['NAXIS1'] = cropped_data.shape[1]
            header['NAXIS2'] = cropped_data.shape[0]

        fits.writeto(output_path, cropped_data, header, overwrite=True)
        print(f"\nФайл успешно обрезан и сохранён как {output_path}")
        print(f"Новый размер: {cropped_data.shape}")


def estimate_source_radius(x, y, threshold, data, max_radius=5):
    x, y = int(round(x)), int(round(y))
    center_value = data[y, x]
    neighbors = [
        data[y+1, x], data[y-1, x],  # верхний и нижний
        data[y, x+1], data[y, x-1]    # правый и левый
    ]
    print(threshold, center_value)
    for i in neighbors:
      if 4*threshold > i:
        return 1

    radii = []


    for dx, dy in [(1,0), (-1,0), (0,1), (0,-1)]:
        r = 1
        while r <= max_radius:
            xi, yi = x + dx*r, y + dy*r
            if xi < 0 or xi >= data.shape[1] or yi < 0 or yi >= data.shape[0]:
                break

            if data[yi, xi] <= center_value*0.5:
                radii.append(r)
                break

            r += 2
        else:
            radii.append(max_radius)
    if np.isfinite(np.mean(radii)):
        return np.mean(radii) if radii else 1
    else:
        return 1



def visualize_source(fits_file, sources):
    with fits.open(fits_file) as hdul:
        data = hdul[0].data

    plt.figure(figsize=(10, 10))
    plt.imshow(data, cmap='gray', norm=LogNorm())

    for x, y, r in sources:
        circle = plt.Circle((x, y), r, color='r', fill=False, linewidth=1)
        plt.gca().add_patch(circle)
        plt.plot(x, y, 'rx', markersize=5)

    plt.title(f'Найдено {len(sources)} источников')
    plt.colorbar(label='Яркость')
    plt.show()


def visualize_peak_profile(x, y, smoothed_data):
    plt.figure(figsize=(10, 4))
    plt.subplot(121)
    plt.plot(np.arange(x - 5, x + 6), smoothed_data[y, x - 5:x + 6], 'b-')
    plt.title(f'X профиль (y={y})')
    plt.xlabel('X координата')
    plt.ylabel('Поток')
    plt.subplot(122)
    plt.plot(np.arange(y - 5, y + 6), smoothed_data[y - 5:y + 6, x], 'g-')
    plt.title(f'Y профиль (x={x})')
    plt.xlabel('Y координата')
    plt.ylabel('Поток')

    plt.tight_layout()
    plt.show()



def find_and_group_sources(fits_file, max_radius=20, ns=0.01):
    with fits.open(fits_file) as hdul:
        data = hdul[0].data

    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    threshold = median + ns *std

    print(f"Статистика фона: median={median:.2f}, std={std:.2f}")
    print(f"Порог обнаружения: {threshold:.2f} (>{ns}σ над фоном)")
    smoothed = gaussian_filter(data, sigma=1)

    peaks = []
    for y in range(5, data.shape[0]-5):
        for x in range(5, data.shape[1]-5):
            if smoothed[y, x] == smoothed[y-5:y+6, x-5:x+6].max() and smoothed[y, x] > threshold:
                region = smoothed[y-5:y+6, x-5:x+6]
                edge_mask = np.zeros((11, 11), dtype=bool)
                edge_mask[1:10, 1:10] = True
                edge_flux = np.mean(region[~edge_mask])
                center_flux = smoothed[y, x]
                edge_mean = np.mean(edge_flux)
                center_flux = smoothed[y, x]
                if edge_flux < 4 * threshold and center_flux > 1.4 * edge_mean:
                    radius = estimate_source_radius(x, y, threshold, smoothed)
                    peaks.append([x, y, radius])
    peaks = np.array(peaks)
    return [(x, y, z) for x, y, z in peaks]

def visualize_results(fits_file, sources):
    with fits.open(fits_file) as hdul:
        data = hdul[0].data

    plt.figure(figsize=(10, 10))
    plt.imshow(data, cmap='gray', norm=LogNorm())

    for x, y, r in sources:
        circle = plt.Circle((x, y), r, color='r', fill=False, linewidth=1)
        plt.gca().add_patch(circle)
        plt.plot(x, y, 'rx', markersize=5)

    plt.title(f'Найдено {len(sources)} источников')
    plt.colorbar(label='Яркость')
    plt.show()


def visualize_source_profile(data, x, y, fwhm, max_analysis_radius, source_num):
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.imshow(data[int(y) - max_analysis_radius:int(y) + max_analysis_radius,
               int(x) - max_analysis_radius:int(x) + max_analysis_radius],
               cmap='gray')
    plt.scatter(max_analysis_radius, max_analysis_radius, c='red', marker='+')
    plt.title(f'Источник {source_num}: ({x:.1f}, {y:.1f})')

    yy, xx = np.indices(data.shape)
    distances = np.sqrt((xx - x) ** 2 + (yy - y) ** 2).astype(int)

    radial_profile = []
    for r in range(0, max_analysis_radius + 1):
        mask = (distances == r)
        if np.any(mask):
            mean_intensity = np.mean(data[mask])
            radial_profile.append((r, mean_intensity))

    if radial_profile:
        radii, intensities = zip(*radial_profile)
        intensities_norm = (np.array(intensities) - np.min(intensities)) / (np.max(intensities) - np.min(intensities))

        plt.subplot(1, 2, 2)
        plt.plot(radii, intensities_norm, 'bo-')
        plt.axhline(y=0.5, color='r', linestyle='--')
        plt.axvline(x=fwhm / 2, color='g', linestyle=':')
        plt.xlabel('Радиус (пиксели)')
        plt.ylabel('Нормализованная интенсивность')
        plt.title(f'Профиль (FWHM = {fwhm:.2f} px)')
        plt.grid(True)

    plt.tight_layout()
    plt.show()


def calculate_fwhm_for_sources(fits_file: str, sources: list, max_analysis_radius: int = 20, show_plots=False):
    with fits.open(fits_file) as hdul:
        data = hdul[0].data.astype(float)

    results = []

    for i, (x, y, r) in enumerate(sources, 1):
        y_idx, x_idx = int(round(y)), int(round(x))
        yy, xx = np.indices(data.shape)
        distances = np.sqrt((xx - x_idx) ** 2 + (yy - y_idx) ** 2)
        distances = distances.astype(int)

        radial_profile = []
        for radius in range(0, max_analysis_radius + 1):
            mask = (distances == radius)
            if np.any(mask):
                mean_intensity = np.mean(data[mask])
                radial_profile.append((radius, mean_intensity))

        if not radial_profile:
            print(f"Для источника {i} ({x:.1f}, {y:.1f}) не удалось построить профиль!")
            results.append((x, y, r, np.nan))
            continue

        radii, intensities = zip(*radial_profile)
        radii = np.array(radii)
        intensities = np.array(intensities)
        intensities_norm = (intensities - np.min(intensities)) / (np.max(intensities) - np.min(intensities))

        try:
            half_max = 0.5
            f = interp1d(radii, intensities_norm - half_max, kind='linear')
            r_half = brentq(f, 0, max_analysis_radius)
            fwhm = 2 * r_half
        except Exception as e:
            print(f"Ошибка вычисления FWHM для источника {i}: {str(e)}")
            fwhm = 2

        results.append((x, y, r, fwhm))

        if show_plots:
            visualize_source_profile(data, x, y, r, fwhm, max_analysis_radius, i)

    return results


def subtract_objects_with_psf(input_file, sources, output_file=None, show_plots=True,
                            subtraction_factor=0.8):

    with fits.open(input_file) as hdul:
        data = hdul[0].data.astype(float)
        header = hdul[0].header

    original_data = data.copy()
    cleaned_data = data.copy()
    psf_model_total = np.zeros_like(data)
    def moffat_2d(A, a, b, c, d, x0, y0):
        return lambda y, x: A*(1 + a*(x-x0)**2 + b*(y-y0)**2 + c*(x-x0)*(y-y0))**-d
    if show_plots:
        plt.figure(figsize=(15, 5))
        plt.ion()  # Включение интерактивного режима

    for i, (x, y, r, fwhm) in enumerate(sources):
        radius = int(np.ceil(fwhm))
        x, y = int(round(x)), int(round(y))

        if (y < radius or y >= data.shape[0]-radius or
            x < radius or x >= data.shape[1]-radius):
            continue

        roi = cleaned_data[y-radius:y+radius, x-radius:x+radius]
        if show_plots:
            plt.clf()
            plt.subplot(1, 3, 1)
            plt.imshow(original_data, cmap='gray',
                      vmin=np.percentile(original_data, 5),
                      vmax=np.percentile(original_data, 95))
            circle = plt.Circle((x, y), radius, color='r', fill=False, linewidth=2)
            plt.gca().add_patch(circle)
            plt.title(f'Оригинал (звезда {i+1}/{len(sources)})')

        initial_params = (roi.max(), 0.1, 0.1, 0.0, 2.5, radius, radius)

        def error_func(p):
            model = moffat_2d(*p)(*np.indices(roi.shape))
            return np.ravel(model - roi)

        fitted_params, _ = leastsq(error_func, initial_params, maxfev=1000, ftol=0.05)
        psf = moffat_2d(*fitted_params)(*np.indices(roi.shape))

        if show_plots:
            plt.subplot(1, 3, 2)
            plt.imshow(psf, cmap='gray')
            plt.title(f'Модель PSF (FWHM={fwhm:.1f})')

        cleaned_roi = roi - psf * subtraction_factor
        background = np.median(roi[roi < np.percentile(roi, 70)])
        cleaned_roi = np.clip(cleaned_roi, background, 0.7*np.mean(roi))

        cleaned_data[y-radius:y+radius, x-radius:x+radius] = cleaned_roi
        psf_model_total[y-radius:y+radius, x-radius:x+radius] += psf
        if show_plots:
            plt.subplot(1, 3, 3)
            plt.imshow(cleaned_data, cmap='gray',
                      vmin=np.percentile(original_data, 5),
                      vmax=np.percentile(original_data, 95))
            plt.title(f'После вычитания {i+1} звезд')
            plt.draw()

    if show_plots:
        plt.ioff()
        plt.figure(figsize=(18, 6))

        plt.subplot(1, 3, 1)
        plt.imshow(original_data, cmap='gray',
                  vmin=np.percentile(original_data, 5),
                  vmax=np.percentile(original_data, 95))
        plt.title('Оригинальное изображение')

        plt.subplot(1, 3, 2)
        plt.imshow(psf_model_total, cmap='gray')
        plt.title('Совокупная модель PSF')

        plt.subplot(1, 3, 3)
        plt.imshow(cleaned_data, cmap='gray',
                  vmin=np.percentile(original_data, 5),
                  vmax=np.percentile(original_data, 95))
        plt.title('Финальный результат')

        plt.tight_layout()
        plt.show()

    if output_file:
        fits.writeto(output_file, cleaned_data, header, overwrite=True)

    return cleaned_data, psf_model_total



def subtract_objects_with_psf(input_file, sources, output_file=None, show_plots=True,
                            subtraction_factor=0.8, smoothing_sigma=1.5):
    """
    Версия с визуализацией на полном изображении и отметкой текущей звезды.
    """
    # Загрузка данных
    with fits.open(input_file) as hdul:
        data = hdul[0].data.astype(float)
        header = hdul[0].header

    original_data = data.copy()
    cleaned_data = data.copy()
    psf_model_total = np.zeros_like(data)

    # Модифицированная функция Moffat
    def moffat_2d(A, a, b, c, d, x0, y0):
        return lambda y, x: A*(1 + a*(x-x0)**2 + b*(y-y0)**2 + c*(x-x0)*(y-y0))**-d

    # Создаем фигуру для анимации процесса
    if show_plots:
        plt.figure(figsize=(15, 5))
        plt.ion()  # Включение интерактивного режима

    for i, (x, y, r, fwhm) in enumerate(sources):
        radius = int(np.ceil(fwhm))
        x, y = int(round(x)), int(round(y))

        if (y < radius or y >= data.shape[0]-radius or
            x < radius or x >= data.shape[1]-radius):
            continue

        roi = cleaned_data[y-radius:y+radius, x-radius:x+radius]

        # Визуализация ДО обработки на полном изображении
        if show_plots:
            plt.clf()
            plt.subplot(1, 3, 1)
            plt.imshow(original_data, cmap='gray',
                      vmin=np.percentile(original_data, 5),
                      vmax=np.percentile(original_data, 95))
            circle = plt.Circle((x, y), radius, color='r', fill=False, linewidth=2)
            plt.gca().add_patch(circle)
            plt.title(f'Оригинал (звезда {i+1}/{len(sources)})')

        initial_params = (roi.max(), 0.1, 0.1, 0.0, 2.5, radius, radius)

        def error_func(p):
            model = moffat_2d(*p)(*np.indices(roi.shape))
            return np.ravel(model - roi)

        fitted_params, _ = leastsq(error_func, initial_params, maxfev=1000, ftol=0.05)
        psf = moffat_2d(*fitted_params)(*np.indices(roi.shape))

        # Визуализация модели PSF
        if show_plots:
            plt.subplot(1, 3, 2)
            plt.imshow(psf, cmap='gray')
            plt.title(f'Модель PSF (FWHM={fwhm:.1f})')

        # Вычитание и обработка
        cleaned_roi = roi - psf * subtraction_factor
        background = np.median(roi[roi < np.percentile(roi, 70)])
        cleaned_roi = np.clip(cleaned_roi, background, 0.7*np.mean(roi))

        cleaned_data[y-radius:y+radius, x-radius:x+radius] = cleaned_roi
        psf_model_total[y-radius:y+radius, x-radius:x+radius] += psf

        # Визуализация ПОСЛЕ обработки
        if show_plots:
            plt.subplot(1, 3, 3)
            plt.imshow(cleaned_data, cmap='gray',
                      vmin=np.percentile(original_data, 5),
                      vmax=np.percentile(original_data, 95))
            plt.title(f'После вычитания {i+1} звезд')
            plt.draw()
            plt.pause(0.5)  # Пауза для просмотра

    # Финальная визуализация
    if show_plots:
        plt.ioff()  # Выключение интерактивного режима
        plt.figure(figsize=(18, 6))

        plt.subplot(1, 3, 1)
        plt.imshow(original_data, cmap='gray',
                  vmin=np.percentile(original_data, 5),
                  vmax=np.percentile(original_data, 95))
        plt.title('Оригинальное изображение')

        plt.subplot(1, 3, 2)
        plt.imshow(psf_model_total, cmap='gray')
        plt.title('Совокупная модель PSF')

        plt.subplot(1, 3, 3)
        plt.imshow(cleaned_data, cmap='gray',
                  vmin=np.percentile(original_data, 5),
                  vmax=np.percentile(original_data, 95))
        plt.title('Финальный результат')

        plt.tight_layout()
        plt.show()

    if output_file:
        fits.writeto(output_file, cleaned_data, header, overwrite=True)

    return cleaned_data, psf_model_total