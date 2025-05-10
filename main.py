from autofiz.commands import visualize, crop, find_sources, subtract
import os


def main():
    current_file = None
    sources = None

    def load_file():
        while True:
            file_path = input("Введите путь к FITS-файлу для анализа: ").strip()
            if os.path.exists(file_path):
                return file_path
            print("Ошибка: файл не найден! Попробуйте снова.")

    current_file = load_file()
    print("\nФайл успешно загружен!")

    while True:
        print("\nТекущий файл:", current_file)


        print("\nДоступные команды:")
        print("1. Визуализировать текущий FITS-файл")
        print("2. Обрезать текущий FITS-файл")
        print("3. Найти точечные объекты в файле")
        print("4. Вычесть найденные объекты из изображения")
        print("5. Сменить анализируемый файл")
        print("0. Выход")

        command = input("\nВведите номер команды (0-5): ").strip()

        if command == '0':
            print("Выход из программы.")
            break

        elif command == '1':
            visualize.run(current_file)

        elif command == '2':
            current_file = crop.run(current_file)
            sources = None  # Сбрасываем найденные источники после обрезки

        elif command == '3':
            sources = find_sources.run(current_file)

        elif command == '4':
            if sources is None:
                sources = find_sources.run(current_file)
            subtract.run(current_file, sources)

        elif command == '5':
            new_file = load_file()
            if new_file != current_file:
                current_file = new_file
                sources = None  # Сбрасываем источники при смене файла
                print("\nФайл успешно изменён!")

        else:
            print("Неизвестная команда. Попробуйте снова.")


if __name__ == "__main__":
    main()
