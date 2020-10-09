#!/bin/sh

GOOGLETEST_VERSION=1.10.0

# googletest

rm -rf $BASEDIR/googletest*
curl -L https://github.com/google/googletest/archive/release-${GOOGLETEST_VERSION}.tar.gz | tar zxf - -C ${BASEDIR}
mv -f ${BASEDIR}/googletest-release-${GOOGLETEST_VERSION}/googletest ${BASEDIR}
rm -rf ${BASEDIR}/googletest/docs ${BASEDIR}/googletest/samples ${BASEDIR}/googletest/scripts ${BASEDIR}/googletest/test 
rm -rf ${BASEDIR}/googletest-release-${GOOGLETEST_VERSION}
