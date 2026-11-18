"""A stand-in for MYSTRAN (test_solve.py, test_connector.py): it plays a run from the deck's first
line, ``$ fake: ok|warn|fatal|fatal989|empty|crash|sleep``, writing the files MYSTRAN 19 writes
(upper-case names, a 144-byte OP2 on a fatal, exit code 0 either way). No pytest import: it runs as
a script, started by the launcher :func:`fake_solver` writes.
"""
import sys
import time
from pathlib import Path

BEGIN = "\n MYSTRAN Version 19.0.0\n\n >> MYSTRAN BEGIN  : 11/17/2026 at 19:50: 1.373 The input file is {name}\n\n"
END = "\n >> MYSTRAN END    : 11/17/2026 at 19:50:38.940\n\n                 * * * END OF JOB * * *\n"
MESSAGES = {
    "ok": " *INFORMATION: NUMBER OF GRID POINTS                                    =         1025\n",
    "warn": (" *WARNING    : GRID        7 HAS NO ELEMENTS ATTACHED                      \n"
             "               IT IS CONSTRAINED BY AUTOSPC\n"
             "\n"
             " *WARNING    : GRID        7 HAS NO ELEMENTS ATTACHED                      \n"
             "               IT IS CONSTRAINED BY AUTOSPC\n"),
    "fatal": (" *ERROR  1964: PSOLID ENTRY HAD FIELD FOR IN BLANK WITH ELEMENT TYPE HEXA8   \n"
              " Stopping processing due to        1 errors found in subr EMG\n"),
    "fatal989": (" *ERROR   989: THE FACTORIZATION OF THE MATRIX KLL      COULD NOT BE COMPLETED BY LAPACK SUBR DPBSTF  \n"
                 "               THE LEADING MINOR OF ORDER         2981 IS NOT POSITIVE DEFINITE\n"
                 "               THIS IS FOR ROW AND COL IN THE MATRIX FOR GRID POINT     1019 COMPONENT   2\n"),
    "empty": " *INFORMATION: NOTHING WAS REQUESTED\n",
}
# the OP2 MYSTRAN 19.0.0 wrote after a fatal: a header (NASTRAN FORT TAPE ID CODE), no tables
EMPTY_OP2_FILE = bytes.fromhex(
    "0400000003000000040000000c0000000300000018000000150000000c0000000400000007000000040000001c000000"
    "4e41535452414e20464f5254205441504520494420434f4445202d201c00000004000000020000000400000008000000"
    "58585858585858580800000004000000ffffffff04000000040000000000000004000000040000000000000004000000")
EMPTY_OP2 = len(EMPTY_OP2_FILE)  # 144 bytes
SLEEP = 8.0  # seconds the `sleep` run takes: longer than any time limit a test gives it


def main(argv) -> int:
    deck = Path(argv[1])
    mode = deck.read_text().splitlines()[0].partition("fake:")[2].strip()
    if mode == "crash":
        print("forrtl: severe (157): Program Exception - access violation")
        return 3
    if mode == "sleep":
        time.sleep(SLEEP)
        mode = "ok"
    text = BEGIN.format(name=deck.name) + MESSAGES[mode] + END
    deck.with_name(deck.stem + ".F06").write_text(text)
    deck.with_name(deck.stem + ".ERR").write_text(text)
    deck.with_name(deck.stem + ".OP2").write_bytes(b"\0" * 4096 if mode in ("ok", "warn") else EMPTY_OP2_FILE)
    return 0  # MYSTRAN exits 0 on a fatal too


def fake_solver(folder) -> Path:
    """An executable in ``folder`` that runs this script with the running interpreter, as a child
    of the launcher (no ``exec``): a time limit has to stop both, like a solver behind a wrapper."""
    folder = Path(folder)
    if sys.platform.startswith("win"):
        path = folder / "fake_mystran.cmd"
        path.write_text(f'@"{sys.executable}" "{Path(__file__).resolve()}" %*\n')
    else:
        path = folder / "fake_mystran"
        path.write_text(f'#!/bin/sh\n"{sys.executable}" "{Path(__file__).resolve()}" "$@"\n')
        path.chmod(0o755)
    return path


def write_deck(path, mode: str) -> Path:
    """A deck that makes the fake solver play ``mode``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"$ fake: {mode}\nSOL 101\nCEND\nBEGIN BULK\nENDDATA\n")
    return path


if __name__ == "__main__":
    sys.exit(main(sys.argv))
