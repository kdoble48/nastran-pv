"""The ``nastran-pv`` command and ``python -m nastran_pv``.

``--version`` is answered here, before :mod:`nastran_pv.cli` is imported: the CLI imports numpy,
VTK and pyNastran at module level, and the version report must still print when one of them is
missing or broken (:mod:`nastran_pv.versions`). Everything else goes to ``cli.main``.
"""
import sys


def main(argv=None) -> int:
    args = sys.argv[1:] if argv is None else list(argv)
    if args[:1] == ["--version"]:
        from .versions import report

        print(report())
        return 0
    from .cli import main as cli_main

    return cli_main(argv)


if __name__ == "__main__":
    sys.exit(main())
