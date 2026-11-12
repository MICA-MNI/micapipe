"""Fail unless micapipe recorded a module as COMPLETED (or wrote its report).

Called by every rule after its micapipe call; see the Snakefile.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))

import micapipe_workflow as mw  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", required=True)
    parser.add_argument("--module", default="", help="expected Module; empty for a QC report")
    parser.add_argument("--subject", required=True)
    parser.add_argument("--session", required=True)
    parser.add_argument("--log", required=True)
    args = parser.parse_args()

    try:
        if args.module:
            mw.check_module_record(args.record, args.module, args.subject, args.session)
        elif not os.path.isfile(args.record) or os.path.getsize(args.record) == 0:
            raise RuntimeError(f"micapipe wrote no report at {args.record}")
    except RuntimeError as err:
        message = f"[micapipe_snakebids] {err}"
        with open(args.log, "a") as stream:
            stream.write("\n" + message + "\n")
        print(f"{message} (log: {args.log})", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
