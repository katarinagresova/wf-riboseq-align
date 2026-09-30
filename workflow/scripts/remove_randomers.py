"""Strip the 4 nt randomer from each end of every read, appending both to its name.

Port of remove8N_stdinout.pl (Emanuel Wyler; stdin/stdout by Gabriel Villamil).
@name -> @name_<5' randomer>:<3' randomer>.

Usage: remove_randomers.py < in.fastq > out.fastq
"""
import sys


def main():
    lines = (line.rstrip(b"\n") for line in sys.stdin.buffer)
    out = sys.stdout.buffer
    for name, seq, plus, qual in zip(lines, lines, lines, lines):
        out.write(b"%s_%s:%s\n%s\n%s\n%s\n"
                  % (name, seq[:4], seq[-4:], seq[4:-4], plus, qual[4:-4]))


if __name__ == "__main__":
    main()
