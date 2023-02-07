import os


def make_dir(path: str) -> None:
    """
    Makes a directory if not present already

    :param path: Path of the directory
    """
    if not os.path.isdir(path):
        os.mkdir(path)
        name = path.split('/')
        print(f"{name[-1]} created successfully")
