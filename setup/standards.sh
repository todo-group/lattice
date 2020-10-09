#!/bin/sh

STANDARDS_VERSION=0.1

# standards

rm -rf $BASEDIR/standards
curl -L https://github.com/todo-group/standards/archive/${STANDARDS_VERSION}.tar.gz | tar zxf - -C ${BASEDIR}
mv -f ${BASEDIR}/standards-${STANDARDS_VERSION}/standards ${BASEDIR}
mv -f ${BASEDIR}/standards-${STANDARDS_VERSION}/config/* ${BASEDIR}/config
rm -rf ${BASEDIR}/standards-${STANDARDS_VERSION}
