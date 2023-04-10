#!/bin/bash

starch=$1
prefix=$2
SPOT=$3
if [[ -z "$starch" || "$(unstarch --is-starch "$starch")" == "0" ]]; then
  echo "Couldn't read starch file '$starch'" >&2
  exit 1
fi

echo -e "$prefix-num-bases\t$(unstarch --bases    "$starch")"
echo -e "$prefix-num-spots\t$(unstarch --elements "$starch")"

if [[ -n "$SPOT" ]]; then
  if [[ -f "$SPOT" ]]; then
    tmp=$(awk 'END{print $NF}' "$SPOT")
    # Check to make sure we read a number (must accept 0 <= x <= 1)
    if [[ ! "$tmp" =~ ^[0-1]?(.[0-9]+)?$ ]]; then
      echo "Couldn't read SPOT file $SPOT" >&2
      exit 1
    fi
    SPOT=$tmp
  fi
  echo -e "$prefix-SPOT\t$SPOT"
fi