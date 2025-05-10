from autofiz.utils import find_and_group_sources, visualize_source


def run(current_file):
    try:
        max_radius = int(input("Максимальный радиус группы (пикселей, например 10): ") or "10")

        print("\nПоиск источников...")
        sources = find_and_group_sources(current_file, max_radius=max_radius)
        print(f"\nНайдено {len(sources)} источников:")

        for i, (x, y, r) in enumerate(sources, 1):
            print(f"{i}. Центр: ({x:.1f}, {y:.1f}), Радиус: {r:.1f} пикселей")

        show_plot = input("\nПоказать визуализацию источников? (y/n): ").lower()
        if show_plot == 'y':
            visualize_source(current_file, sources)

        return sources

    except Exception as e:
        print(f"Ошибка при поиске источников: {str(e)}")
        return None
