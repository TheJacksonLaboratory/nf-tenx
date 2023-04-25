from argparse import ArgumentParser, Namespace
from pathlib import Path


def parse_cl_paths(*expected_args: tuple[str, str, str]) -> Namespace | None:
    """
    Takes expected command-line arguments (representing paths) and returns
    an `argparse.Namespace` of `pathlib.Path`s.

    Parameters
    ----------
    `*expected_args` : `tuple[str, str, str]`
        `(long_name, short_name, help)`. These are eventually passed into
        `argparse.ArgumentParser.add_argument()` as keyword arguments.
        Note that all command-line arguments will be read as
        `pathlib.Path` objects.

    Returns
    -------
    `args` : `arparse.Namespace`
        The parsed command-line arguments, each accessible by
        `args.long_name`.

    Raises
    ------
    `ValueError`
        If any of the tuples are not length 2 or 3, raise `ValueError`
    """

    # Initialize parser
    parser = ArgumentParser()

    # Do the same for 3-tuples, where final element is help message
    for long_name, short_name, help_msg in (arg for arg in expected_args):
        parser.add_argument(long_name, short_name, help=help_msg, type=Path)

    # Return the resulting Namepsace object
    return parser.parse_args()
