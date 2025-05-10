from autofiz.utils import visualize_fits, intetn_fits


def run(current_file):

    visualize_fits(current_file)

    show_more = input("Показать распределение интенсивностей? (y/n): ").lower()
    if show_more == 'y':
        intetn_fits(current_file)

    return current_file
