RMFLAGS = -f

all: config

config:
	mkdir pickles/
	mkdir reads/
	mkdir samples/
clean:
	rm -rf pickles/
	rm -rf reads
	rm -rf samples
