NIST to MassBank converter
==========================

Compile the converter with

```
mvn clean package
```

run the result with:

```
java -jar target/MassBank2NIST-0.0.1-SNAPSHOT.jar 
```

Usage:

```
[MassBank-to-Library]: -mode mb2lib  -spectra path -mol path -output path
[Library-to-MassBank]: -mode lib2mb  -library file -output path -prefix prefix
```



