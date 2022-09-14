import os


def make_dir(path):
    if not os.path.isdir(path):
        print('hello')
        os.mkdir(path)


if __name__ == '__main__':
    make_dir('Result')
