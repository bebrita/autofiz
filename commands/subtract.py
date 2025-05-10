import os

from autofiz.utils import calculate_fwhm_for_sources, subtract_objects_with_psf, visualize_fits


def run(current_file, sources):
    try:
        print("\nАнализ FWHM для найденных источников...")
        results = calculate_fwhm_for_sources(current_file, sources)

        show_results = input("\nПоказать результаты анализа FWHM? (y/n): ").lower()
        if show_results == 'y':
            for i, (x, y, r, fwhm) in enumerate(results, 1):
                print(f"Источник {i}: X={x:.1f}, Y={y:.1f}, R={r:.1f} → FWHM={fwhm:.2f} пикселей")

        base_name, ext = os.path.splitext(current_file)
        output_path = f"{base_name}_clean{ext}"

        if output_path:
            print("\nВычитание объектов...")
            cleaned, psf_model = subtract_objects_with_psf(
                current_file,
                results,
                output_file=output_path,
                show_plots=False
            )

            print(f"Очищенное изображение сохранено как {output_path}")

            show_final = input("Показать сравнение до/после? (y/n): ").lower()
            if show_final == 'y':
                visualize_fits(current_file)
                visualize_fits(output_path)
        else:
            print("Сохранение отменено.")

    except Exception as e:
        print(f"Ошибка при вычитании объектов: {str(e)}")
