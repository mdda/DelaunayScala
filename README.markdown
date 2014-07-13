# Delaunay Triangulation in Scala

## Motivation

After looking around for a decent stand-alone Scala implementation of Delaunay triangulation, 
I gave up finding anything that seemed idiomatic to me.

The approach is per the "bourke method" described in http://paulbourke.net/papers/triangulate/. 
However, it's heavily modified here, because the original code is very mutable-state heavy, 
whereas this scala version is 'pure immutable'.


## File Layout

This works in a straight-forward ```sbt``` setting.  
The meat of it is just the ```src/main/scala/Delaunay.scala``` file.  

Please let me know if there are insufficient comments.

 
## Usage

```scala
import Deluany
val triangles = Deluanay.Triangulation(measurements.toList)
```


## Dependencies

No dependencies : This is really a single scala file.

One can easily redefine 'Vector2' to be the same as the corresponding class in your own projects : 
It is just a ```Seq[ (Float, Float) ]```.

 
## Installation

The scala file also contains a ```main```, so that you can run a simple demo test.  

```bash
cd DelaunayScala
sbt
run 
```

and then copy-paste the output into a 'Processing' session (free), to see the geometry of the output.

