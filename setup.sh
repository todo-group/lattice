#!/bin/sh

LIBS="googletest standards"

BASEDIR=$(cd "$(dirname $0)"; pwd)
for i in ${LIBS}; do
  echo "[${i}]"
  env BASEDIR=${BASEDIR} sh ${BASEDIR}/config/${i}.sh
done
