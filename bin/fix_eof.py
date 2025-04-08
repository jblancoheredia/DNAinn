import os
import sys


def fix_bam(filename):
    header = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
    eof = (
        "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC"
        "\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    )
    if not os.path.isfile(filename):
        sys.exit("Missing file %s" % filename)
    size = os.path.getsize(filename)
    h = open(filename, "rb")  # read only for now
    # Check it looks like a BGZF file
    # (could still be GZIP'd, in which case the extra block is harmless)
    data = h.read(len(header))
    if data != header:
        sys.exit("File %s is not a BAM file" % filename)
    # Check if it has the EOF already
    h.seek(size - 28)
    data = h.read(28)
    h.close()
    if data == eof:
        sys.stderr.write("EOF already present in %s\n" % filename)
    else:
        sys.stderr.write("Adding EOF block to %s\n" % filename)
        h = open(filename, "ab")
        h.write(eof)
        h.close()


if len(sys.argv) == 1:
    sys.exit("Takes one or more BGZF/BAM filenames as arguments (edits in place)")
for bam_filename in sys.argv[1:]:
    fix_bam(bam_filename)