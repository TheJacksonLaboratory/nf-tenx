from argparse import ArgumentParser, Namespace
from pathlib import Path


def parse_cl_paths(expected_args: dict[str, str]) -> Namespace:
    """
    Takes expected command-line arguments (representing paths) and
    returns an argparse.Namespace of pathlib.Paths.

    Parameters
    ----------
    expected_args : dict[str, str]
        (long_name, short_name). These are eventually passed into
        argparse.ArgumentParser.add_argument() as keyword arguments.
        Note that all command-line arguments will be read as
        pathlib.Path objects.

    Returns
    -------
    args : arparse.Namespace
        The parsed command-line arguments, each accessible by
        args.long_name.
    """

    # Initialize parser
    parser = ArgumentParser()

    # Iterate over each argument's key-value pair, adding to arg_parser
    for long_name, short_name in expected_args.items():
        parser.add_argument(long_name, short_name, type=Path)

    # Return the resulting Namepsace object
    return parser.parse_args()
