from utils.conf import MUTABIND2_URL
from utils.driver_conf import initialize_driver


def main():
    open_default_page()




def open_default_page():
    driver = initialize_driver()
    driver.get(MUTABIND2_URL)


if __name__ == '__main__':
    main()