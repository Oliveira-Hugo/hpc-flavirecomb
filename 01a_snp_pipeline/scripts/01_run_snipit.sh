#!/usr/bin/env bash

set -euo pipefail

FASTA="$1"
OUTDIR="$2"

echo "=== DEBUG ==="
echo "FASTA = ${FASTA}"
echo "OUTDIR = ${OUTDIR}"
echo "HOSTNAME = $(hostname)"
echo "PWD(before cd) = $(pwd)"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "PWD(after cd) = $(pwd)"

snipit -s "${FASTA}"

