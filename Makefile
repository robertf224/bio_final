RMFLAGS = -f

all: config

config:
	mkdir pickles/
	mkdir reads/
	mkdir samples/
	mkdir plots/
clean:
	rm -rf pickles/
	rm -rf reads/
	rm -rf samples/
	rm -rf plots/
