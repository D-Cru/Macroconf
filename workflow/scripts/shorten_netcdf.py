#!/usr/bin/env python
# Taken & adapted from https://gist.github.com/opie4624/3896526

# import modules used here -- sys is a very standard one
import argparse, logging
import mdtraj as md


# Gather our code in a main() function
def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    print("Processing trajectory.")
    logging.info("You passed an argument.")
    logging.debug(f"Topology: {args.topology}")
    logging.debug(f"Input trajectory: {args.input_traj}")
    logging.debug(f"Output trajectory: {args.output_traj}")
    logging.debug(f"Stride: {args.stride}")
    t = md.load(args.input_traj, top=args.topology)
    t[:: args.stride].save_netcdf(args.output_traj)


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Does a thing to some stuff.",
        epilog="As an alternative to the commandline, params can be placed in a file, one per line, and specified on the commandline like '%(prog)s @params.conf'.",
        fromfile_prefix_chars="@",
    )
    # TODO Specify your real parameters here.
    parser.add_argument(
        "topology", help="pass topology to the program", metavar="topology"
    )
    parser.add_argument(
        "input_traj",
        help="pass input .netcdf trajectory to the program",
        metavar="traj-in",
    )
    parser.add_argument(
        "stride",
        help="pass stride for shortening traj to the program",
        metavar="stride",
        type=int,
    )
    parser.add_argument(
        "output_traj",
        help="pass ouput .netcdf trajectory to the program",
        metavar="traj-out",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="increase output verbosity",
        action="store_true",
    )
    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    main(args, loglevel)
