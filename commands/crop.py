from astropy.io import fits
import os
from prog.utils import crop_fits


def run(current_file):

    print(f"Текущий файл: {current_file}")
    print(f"Размер изображения: {fits.getdata(current_file).shape}")

    try:
        data = fits.getdata(current_file)
        height, width = data.shape
        while True:
            x_start = int(input(f"Начало по X (0-{width - 1}): "))
            if 0 <= x_start < width:
                break
            print(f"Ошибка: значение должно быть от 0 до {width - 1}")

        while True:
            x_end = int(input(f"Конец по X ({x_start + 1}-{width}): "))
            if x_start < x_end <= width:
                break
            print(f"Ошибка: значение должно быть от {x_start + 1} до {width}")

        while True:
            y_start = int(input(f"Начало по Y (0-{height - 1}): "))
            if 0 <= y_start < height:
                break
            print(f"Ошибка: значение должно быть от 0 до {height - 1}")

        while True:
            y_end = int(input(f"Конец по Y ({y_start + 1}-{height}): "))
            if y_start < y_end <= height:
                break
            print(f"Ошибка: значение должно быть от {y_start + 1} до {height}")

        base_name, ext = os.path.splitext(current_file)
        output_path = f"{base_name}_cut{ext}"
        crop_fits(current_file, output_path, (x_start, x_end), (y_start, y_end))
        return output_path

    except Exception as e:
        print(f"Ошибка при обрезке: {str(e)}")
        return current_file