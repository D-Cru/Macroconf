"""
Common code for unit testing of rules generated with Snakemake 7.12.1.
"""

from pathlib import Path
import subprocess as sp
import os
import sys
import warnings


class OutputChecker:
    def __init__(self, data_path, expected_path, workdir):
        self.data_path = data_path
        self.expected_path = expected_path
        self.workdir = workdir

    def check(self):
        input_files = set(
            (Path(path) / f).relative_to(self.data_path)
            for path, subdirs, files in os.walk(self.data_path)
            for f in files
        )
        expected_files = set(
            (Path(path) / f).relative_to(self.expected_path)
            for path, subdirs, files in os.walk(self.expected_path)
            for f in files
        )
        unexpected_files = set()
        for path, subdirs, files in os.walk(self.workdir):
            for f in files:
                f = (Path(path) / f).relative_to(self.workdir)
                if str(f).startswith(".snakemake"):
                    continue
                if f in expected_files:
                    self.compare_files(
                        self.workdir / f, self.expected_path / f
                    )
                elif f in input_files:
                    # ignore input files
                    pass
                else:
                    unexpected_files.add(f)
            # # check for missing files
            # missing_files = (expected_files | input_files) & set(files)
        if unexpected_files:
            for fi in unexpected_files:
                if str(fi).endswith(".log") or str(fi).endswith(".ipynb"):
                    warnings.warn("Unexpected file (log-file): {}".format(fi))
                else:
                    raise ValueError(
                        "Unexpected files:\n{}".format(
                            "\n".join(sorted(map(str, unexpected_files)))
                        )
                    )
        # if missing_files:
        #     raise ValueError(
        #         "Missing files:\n{}".format(
        #             "\n".join(sorted(map(str, missing_files)))
        #         )
        #     )

    def compare_files(self, generated_file, expected_file):
        sp.check_output(["cmp", generated_file, expected_file])


class OutputChecker_present(OutputChecker):
    def compare_files(self, generated_file, expected_file):
        pass


class OutputChecker_smart(OutputChecker):
    def compare_files(self, generated_file, expected_file):
        byte_files = set([".yaml", ".txt", ".csv", ".tsv", ".log"])

        if generated_file.suffix in byte_files:
            try:
                sp.check_output(["cmp", generated_file, expected_file])
            except sp.CalledProcessError as e:
                print(
                    "Files differ: {} {}".format(
                        generated_file, expected_file
                    ),
                    file=sys.stderr,
                )
                raise e
        else:
            pass
