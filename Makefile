#!/usr/bin/env make

sophie: sophie.c
	gcc -Ofast -Wall -Wextra -Wconversion -std=c99 sophie.c -osophie

clean:
	rm -f sophie