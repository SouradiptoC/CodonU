import os


def make_dir(path):
    if not os.path.isdir(path):
        os.mkdir(path)


def is_file(path):
    return os.path.isfile(path)


if __name__ == '__main__':
    make_dir('Result')
