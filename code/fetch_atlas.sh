#!/bin/bash
# fetch_atlas.sh -- clean, verified download of the PDAC atlas from Zenodo.
#
# Usage (run inside tmux):
#   export TOK='<zenodo share token>'
#   bash fetch_atlas.sh
#
# Deletes any existing copy first. No resume, no -O games -- the previous
# corruption came from mixing rsync with wget -c on a partially moved file.

set -uo pipefail

DEST=/ix/tbruno/anp795/PDAC_Atlas
FNAME=combined_9_datasets_hierarchical_annotation_seurat.rds
URL="https://zenodo.org/records/21795216/files/${FNAME}?download=1"

EXPECTED_MD5=e9361c20b0776e0d66b7068ca4449a36
EXPECTED_SIZE=51999181312

if [ -z "${TOK:-}" ]; then
  echo "ERROR: TOK is not set. Run:  export TOK='<your zenodo token>'" >&2
  exit 1
fi

mkdir -p "$DEST"
cd "$DEST" || exit 1

echo "=== free space ==="
df -h "$DEST" | tail -1

if [ -e "$FNAME" ]; then
  echo "=== removing existing copy ($(stat -c %s "$FNAME") bytes) ==="
  rm -f "$FNAME"
fi
rm -f wget-log

echo "=== downloading (expect ~1.5 h at 8-10 MB/s) ==="
date
wget --progress=dot:giga --tries=20 --waitretry=15 --timeout=60 \
     -O "$FNAME" "${URL}&token=${TOK}"
RC=$?
date

if [ $RC -ne 0 ]; then
  echo "WARNING: wget exited $RC. File may be incomplete."
  echo "Re-run this script to start over, or resume with:"
  echo "  wget -c --progress=dot:giga -O \"$FNAME\" \"${URL}&token=\$TOK\""
fi

echo
echo "=== verifying ==="
ACTUAL_SIZE=$(stat -c %s "$FNAME")
echo "size     expected: $EXPECTED_SIZE"
echo "size       actual: $ACTUAL_SIZE"
if [ "$ACTUAL_SIZE" != "$EXPECTED_SIZE" ]; then
  echo "RESULT: FAIL (wrong size) -- do not use this file"
  exit 1
fi

echo "computing md5 (a few minutes)..."
ACTUAL_MD5=$(md5sum "$FNAME" | awk '{print $1}')
echo "md5      expected: $EXPECTED_MD5"
echo "md5        actual: $ACTUAL_MD5"

if [ "$ACTUAL_MD5" = "$EXPECTED_MD5" ]; then
  echo
  echo "RESULT: PASS -- file is intact and ready to use"
  ls -lh "$FNAME"
  exit 0
else
  echo
  echo "RESULT: FAIL (md5 mismatch) -- delete and re-run this script"
  exit 1
fi
