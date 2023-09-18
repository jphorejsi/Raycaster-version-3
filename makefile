target: raytracer1c.o
	g++ raytracer1c.o -o target

raycast.o: raytracer1c.cc
	g++ -c raytracer1c.cc

clean:
	rm *.o target