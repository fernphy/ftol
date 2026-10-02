#!/bin/bash
# Delete files from the FTOL input data FigShare deposit (19474316).
#
# Usage: bash figshare_delete_files.sh <id>:<name>:<md5> [<id>:<name>:<md5> ...]
#
# Each file is only deleted if its id, name AND md5 all match what FigShare
# currently lists, so a stale or mistyped id can never remove the wrong file.
# Typical use: removing same-named duplicates that retried uploads leave in
# the pending (unpublished) edit of the deposit before publishing. Find ids
# with: curl -H "Authorization: token $FIGSHARE_TOKEN" \
#   https://api.figshare.com/v2/account/articles/19474316/files
# The token is read from .Renviron and never printed. This script never
# publishes anything.
set -euo pipefail

ART=19474316
API=https://api.figshare.com/v2/account/articles/$ART
TOKEN=$(grep -h '^FIGSHARE_TOKEN' /home/jnitta/ftol/.Renviron | head -1 \
  | sed 's/^[^=]*=//; s/["'"'"']//g')

[ $# -gt 0 ] || { echo "usage: $0 id:name:md5 ..." >&2; exit 1; }

list_files() {
  curl -sf -H "Authorization: token $TOKEN" "$API/files" | python3 -c '
import sys, json
for f in json.load(sys.stdin):
    print(f["id"], f["name"], f["size"], f["computed_md5"])'
}

current=$(list_files)

for entry in "$@"; do
  IFS=: read -r id name md5 <<< "$entry"
  if grep -qF "$id $name " <<< "$current" \
     && grep -q "^$id $name .* $md5\$" <<< "$current"; then
    code=$(curl -s -o /dev/null -w '%{http_code}' -X DELETE \
      -H "Authorization: token $TOKEN" "$API/files/$id")
    echo "deleted $id ($name): HTTP $code"
  else
    echo "SKIP $id ($name): not found or name/md5 mismatch"
  fi
done

echo "--- remaining files ---"
list_files
