#!/bin/sh
# Which algebra does each deposit attribute to Maxwell?
#
# This is the one script in the repository that needs the network.  It downloads
# the ten Zenodo deposits by this author, converts each to text with pdftotext,
# and counts the keywords that the attribution claim rests on.  Nothing is
# cached in the repository; rerunning it reproduces the table in PAPER.md Sec. 6
# against whatever is on Zenodo at the time.
#
# Requirements: curl, pdftotext (poppler-utils), python3.
set -e

OUT="${1:-$(mktemp -d)}"
mkdir -p "$OUT"
echo "working directory: $OUT"

# record id : short name          (resolved from the Zenodo REST API, oldest first)
RECORDS="
22876757:acoustic-projector
22877345:x-ternary
22877912:mppt-v1.0
22914238:mppt-v1.1
22923154:cobar
22934058:momentum-drive
22959886:cqft
22962188:pctp
22962557:sotp
23112560:openql-openol
"

for rec in $RECORDS; do
    id=${rec%%:*}; name=${rec##*:}
    if [ ! -f "$OUT/$name.txt" ]; then
        url=$(curl -sS "https://zenodo.org/api/records/$id" \
              | python3 -c 'import json,sys; print(json.load(sys.stdin)["files"][0]["links"]["self"])')
        curl -sS -L "$url" -o "$OUT/$name.pdf"
        pdftotext -layout "$OUT/$name.pdf" "$OUT/$name.txt"
    fi
done

printf '\n%-22s %8s %8s %8s %8s %8s %8s %8s\n' \
       deposit Maxwell octonion quaternion Treatise 1873 1865 Heaviside
printf '%s\n' '----------------------------------------------------------------------------------------------'
for rec in $RECORDS; do
    name=${rec##*:}
    f="$OUT/$name.txt"
    printf '%-22s' "$name"
    for kw in Maxwell octonion quaternion Treatise 1873 1865 Heaviside; do
        printf ' %8d' "$(grep -oiE "$kw" "$f" | wc -l)"
    done
    printf '\n'
done

printf '\nEvery sentence in which Maxwell is named (line breaks folded first):\n'
for rec in $RECORDS; do
    name=${rec##*:}
    tr '\n' ' ' < "$OUT/$name.txt" | tr -s ' ' \
        | grep -oiE '[^.]*Maxwell[^.]*\.' 2>/dev/null \
        | fold -s -w 76 | sed "1s/^/  [$name] /; 2,\$s/^/      /" || true
done
