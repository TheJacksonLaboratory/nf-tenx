from argparse import ArgumentParser, Namespace


def parse_cl(*expected_args: tuple[str, str, type]) -> Namespace:
    """
    Takes expected command-line arguments and returns an
    argparse.Namespace with them.

    Parameters
    ----------
    *expected_args : tuple[str, str, type]
        (long_name, short_name, type). These are eventually passed into
        argparse.ArgumentParser.add_argument() as keyword arguments.

    Returns
    -------
    args : arparse.Namespace
        The parsed command-line arguments, each accessible by
        args.long_name.
    """

    # Initialize parser
    parser = ArgumentParser()

    # Iterate over each argument's key-value pair, adding to arg_parser
    for long_name, short_name, t in expected_args:
        parser.add_argument(long_name, short_name, type=t)

    # Return the resulting Namespace object
    return parser.parse_args()
