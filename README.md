NIST to MassBank converter
==========================

It is required to install the included files into your local maven repository upfront:
```
mvn install:install-file -Dfile=./lib/cdk-1.3.5.jar -DgroupId=org.openscience.cdk  -DartifactId=cdk-bundle -Dversion=1.3.5 -Dpackaging=jar
mvn install:install-file -Dfile=./lib/massbank.jar -DgroupId=massbank  -DartifactId=massbank -Dversion=1.0 -Dpackaging=jar
mvn install:install-file -Dfile=./lib/metfusion.jar -DgroupId=de.ipbhalle.msbi  -DartifactId=metfusion -Dversion=1.0 -Dpackaging=jar
```

Then compile the converter with

```
mvn clean package
```

run the result with:

```
java -jar target/MassBank2NIST-0.0.1-SNAPSHOT.jar 
```

Usage:

```
java -jar target/MassBank2NIST-0.0.1-SNAPSHOT.jar -mode mb2lib  -spectra path -mol path -output path
java -jar target/MassBank2NIST-0.0.1-SNAPSHOT.jar -mode lib2mb  -library primidone.library -output records -prefix SU -properties example.properties
```



