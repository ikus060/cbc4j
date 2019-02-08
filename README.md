[![Build](http://git.patrikdufresne.com/pdsl/cbc4j/badges/master/build.svg)](http://git.patrikdufresne.com/pdsl/cbc4j)

# cbc4j

This project provide a pre-compiled version of [cbc](https://www.coin-or.org/)
for Java. This library is working on linux and windows.

## Compile
This section describe how to configure a linux debian computer to compile the
cbc4j project. A special configuration is required to cross-compile for
windows.

    sudo apt-get install default-jdk maven swig build-essential gcc-mingw-w64 mingw-w64 g++-multilib
    mvn -Drevision=1.0-SNAPSHOT clean install

## Usage
To make use of this library, you must add the following to your maven project.

    <repositories>
        <repository>
            <id>patrikdufresne</id>
            <url>http://nexus.patrikdufresne.com/content/groups/public/</url>
        </repository>
    </repositories>

    ...

    <profile>
        <id>linux_x86_64</id>
        <activation>
            <os>
                <name>linux</name>
                <arch>amd64</arch>
            </os>
        </activation>
        <dependencies>
            <dependency>
                <groupId>com.patrikdufresne.cbc4j</groupId>
                <artifactId>cbc4j-linux-x86_64</artifactId>
                <version>1.1</version>
            </dependency>
        </dependencies>
    </profile>
    <profile>
        <id>win_x86</id>
        <activation>
            <os>
                <family>windows</family>
                <arch>amd64</arch>
            </os>
        </activation>
        <dependencies>
            <dependency>
                <groupId>com.patrikdufresne.cbc4j</groupId>
                <artifactId>cbc4j-win-x86_64</artifactId>
                <version>1.1</version>
            </dependency>
        </dependencies>
    </profile>
