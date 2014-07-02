//package Delaunay

import VSFM.{Vector2,Vector3}

package object Delaunay {
  def Triangulation_n4(measurements : List[Vector2]) : Seq[(Int, Int, Int)] = {
    val n = measurements.length
    
    val x = measurements map { _.x }
    val y = measurements map { _.y }
    val z = (x zip y) map { case (a,b) => (a*a + b*b) }
    
    def is_addable(i:Int,j:Int,k:Int) : Boolean = {
      val xn = (y(j)-y(i)) * (z(k)-z(i)) - (y(k)-y(i)) * (z(j)-z(i))
      val yn = (z(j)-z(i)) * (x(k)-x(i)) - (z(k)-z(i)) * (x(j)-x(i))
      val zn = (x(j)-x(i)) * (y(k)-y(i)) - (x(k)-x(i)) * (y(j)-y(i))

      val flag = (zn<0) && (0 until n).forall { 
         m => (((x(m)-x(i))*xn + (y(m)-y(i))*yn + (z(m)-z(i))*zn) <= 0)
      }
      flag
    }
    
    for {
      i <- (0) until (n-2)
      j <- (i+1) until n
      k <- (i+1) until n
      if (j!=k)
      if is_addable(i,j,k)
    } yield {
      printf(s"Adding (${i}, ${j}, ${k})\n")
      (i,j,k)
    }
  }
  
  /* Do something better than n^4 here... */
  
  /*
    Takes as input NV vertices in array pxyz
    The vertex array pxyz must be big enough to hold 3 more points
    The vertex array must be sorted in increasing x values say
    
    Returned is a list of ntri triangular faces in the array v
    These triangles are arranged in a consistent clockwise order.
    The triangle array 'v' should be malloced to 3 * nv

    
    qsort(p,nv,sizeof(XYZ),XYZCompare);
    
    int XYZCompare(void *v1,void *v2) {
      XYZ *p1,*p2;
      p1 = v1;
      p2 = v2;
      if (p1->x < p2->x)
        return(-1);
      else if (p1->x > p2->x)
        return(1);
      else
        return(0);
    }
  */

  def Triangulation_bourke(measurements : List[Vector2]) : Seq[(Int, Int, Int)] = { // NB: List must be sorted in increasing x...
    case class ITRIANGLE(p1:Int, p2:Int, p3:Int, completed:Boolean)
    case class IEDGE(p1:Int, p2:Int)
    
    class AnihilationSet[IEDGE] extends Set[IEDGE] {
      def add_with_anihilation(that: IEDGE):AnihilationSet[IEDGE] = {
        if( this.contains(that)) {
          this.remove(that)
        }
        else {
          val that_reversed = IEDGE(that.p2, that.p1)
          if( this.contains(that_reversed) ) {
            this.remove(that_reversed)
          }
          else {
            this.add(that)
          }
        }
      }
    }
    
    /*
     {
      override def equals(o: Any) = o match {
        case that: IEDGE(q1, q2) => ((p1==q1) && (p2==q2)) || ((p1==q2) && (p2==q1))
        case _ => false
      }
    }
    */
    
    //case class XYZ(p1:Float, p2:Float, p3:Float)

    val EPSILON = 0.000001

    //  Return TRUE if a point q(x,y) is inside the circumcircle made up of the points p1(x,y), p2(x,y), p3(x,y)
    //  The circumcircle centre (x,y) is returned and the radius r
    //  NOTE: A point on the edge is inside the circumcircle
    def CircumCircle( q:Vector2, p1:Vector2, p2:Vector2, p3:Vector2) : (Boolean, Vector2, Float) {
      if ( Math.abs(p1.y-p2.y) < EPSILON && Math.abs(p2.y-p3.y) < EPSILON ) {
        System.out.println("CircumCircle: Points are colinear");
        (false, new Vector2(0,0), 0)
      }
      else {
        val (mx1, my1) = ( (p1.x+p2.x)/2.0, (p1.y+p2.y)/2.0 )
        val (mx2, my2) = ( (p2.x+p3.x)/2.0, (p2.y+p3.y)/2.0 )
        
        val c = 
          if ( Math.abs(p2.y-p1.y) < EPSILON ) {
            val m2 = -(p3.x-p2.x) / (p3.y-p2.y)
            val xc =  (p2.x+p1.x) / 2.0
            val yc =  m2 * (xc - mx2) + my2
            new Vector2(xc, yc)
          }
          else 
            if ( Math.abs(p3.y-p2.y) < EPSILON ) {
              val m1 = -(p2.x-p1.x) / (p2.y-p1.y)
              val xc =  (p3.x + p2.x) / 2.0
              val yc =  m1 * (xc - mx1) + my1
              new Vector2(xc, yc)
            }
            else {
              val m1 = -(p2.x-p1.x) / (p2.y-p1.y)
              val m2 = -(p3.x-p2.x) / (p3.y-p2.y)
              val xc =  (m1*mx1 - m2*mx2 + my2 - my1) / (m1 - m2)
              val yc =  m1 * (xc - mx1) + my1
              new Vector2(xc, yc)
            }
          
          val rsqr = {
            val (dx, dy) = (p2.x-c.x, p2.y-c.y)
            dx*dx + dy*dy
          }
          val drsqr = {
            val (dx, dy) = (q.x-c.y, q.y-c.y)
            dx*dx + dy*dy
          }
          
          ( drsqr <= rsqr, c, Math.sqrt(rsqr) )
        }
      }

      val n_points = measurements.length // Was nv
      
      // Find the maximum and minimum vertex bounds, to allow calculation of the bounding triangle
      val Pmin = Vector2( measurements.map(_.x).min, measurements.map(_.y).min )  // Top Left
      val Pmax = Vector2( measurements.map(_.x).max, measurements.map(_.y).max )  // Bottom Right
      val diameter = {
        val pdiff = Vector2( Pmax.x - Pmin.x, Pmax.y - Pmin.y)
        (pdiff.x > pdiff.y) ? pdiff.x : pdiff.y
      }
      val Pmid = Vector2( (Pmin.x + Pmax.x)/2.0, (Pmin.y + Pmax.y)/2.0 )
    
      /*
        Set up the supertriangle, which is a triangle which encompasses all the sample points.
        The supertriangle coordinates are added to the end of the vertex list. 
        The supertriangle is the first triangle in the triangle list.
      */
      
      val point = measurements ::: List( 
        Vector2( Pmid.x - 2.0*diameter, Pmid.y - 1.0*diameter), 
        Vector2( Pmid.x - 0.0*diameter, Pmid.y + 2.0*diameter), 
        Vector2( Pmid.x + 2.0*diameter, Pmid.y - 1.0*diameter)
      )
     
      val current_triangles   = List( ITRIANGLE(n_points+0, n_points+1, n_points+2) ) // initially containing the supertriangle
      val completed_triangles = List[ITRIANGLE]                                       // initially empty 
      
      // Add each (original) point, one at a time, into the existing mesh
      for (int i=0; i<n_points; i++) {
        val point_being_added = point(i)
        
        val (completed_triangles_updated : List[ITRIANGLE], current_triangles_updated: List[ITRIANGLE], edges_created: AnihilationSet[IEDGE]) = 
          // Process all the triangles in the current_triangles list
          current_triangles.foldLeft(completed_triangles, Nil, Nil) {
            ((completed, current, edges), triangle) => {
              // If the point 'point_being_added' lies inside the circumcircle then the three edges 
              // of that triangle are added to the edge buffer and that triangle is removed.
              
              // Find the coordinates of the points in this incomplete triangle
              val corner1 = point( triangle.p1 )
              val corner2 = point( triangle.p2 )
              val corner3 = point( triangle.p3 )
              
              val (inside, circle, r) = CircumCircle(point_being_added,  corner1,  corner2,  corner3)
              
              // have we moved too far in x to bother with this one ever again? (initial point list must be sorted for this to work)
              if (circle.x + r < point_being_added.x) {
                ( triangle::completed, current, edges ) // Add this triangle to the 'completed' accumulator, and don't add it on current list
              }
              else {
                if(inside) {
                  // Add the triangle's edge onto the edge pile, and remove the triangle
                  val edges_with_triangle_added = 
                    edges
                      .add_with_anihilation(IEDGE(triangle.p1, triangle,p2))
                      .add_with_anihilation(IEDGE(triangle.p2, triangle,p3))
                      .add_with_anihilation(IEDGE(triangle.p3, triangle,p1))
                  ( completed, current, edges_with_triangle_added )
                }
                else {
                  ( completed, triangle::current, edges )  // This point was not inside this triangle - just add it to the 'current' list
                }
              }
            }
          }
          
      

/*        
        //def add_point_to_triangles(live: 
        
        //val ( live_triangles, completed_triangles, edges ) =  add_point_to_triangles(...)
        
        
        
        // If the point 'point_being_added' lies inside the circumcircle then the three edges 
        // of that triangle are added to the edge buffer and that triangle is removed.
        // Do this for all relevant triangles
        for { 
          triangle <- current_triangles 
          if(!triangle.completed)
        } {
          // Find the coordinates of the points in this incomplete triangle
          val corner1 = point( triangle.p1 )
          val corner2 = point( triangle.p2 )
          val corner3 = point( triangle.p3 )
          
          val (inside, circle, r) = CircumCircle(point_being_added,  corner1,  corner2,  corner3)
          
          // have we moved too far in x to bother with this one ever again? (initial point list must be sorted)
          if (circle.x + r < point_being_added.x) {
            is_triangle_complete(j) = true 
            continue
            // inside is false in this case, for sure
          }
          
          
        }
          if (inside) {
            // Add these edges onto the edge pile
            {
            edges[nedge+0].p1 = v[j].p1;
            edges[nedge+0].p2 = v[j].p2;
            
            edges[nedge+1].p1 = v[j].p2;
            edges[nedge+1].p2 = v[j].p3;
            
            edges[nedge+2].p1 = v[j].p3;
            edges[nedge+2].p2 = v[j].p1;
            
            nedge += 3;
            }
            
            // Delete this triangle (cover it over with another one)
            {
            v[j].p1 = v[ntri-1].p1; 
            v[j].p2 = v[ntri-1].p2;
            v[j].p3 = v[ntri-1].p3;
            
            complete[j] = complete[ntri-1];
            
            ntri--;
            }
            
            // Do this iteration again (one the data we've just swapped in)
            j--;
          }
        }
        
        // Tag multiple edges
        // Note: if all triangles are specified anticlockwise then all
        // interior edges are opposite pointing in direction.
        
        // Essentially, this takes the edge list, and makes duplicates anihilate each other...
        
        for (int j=0;j<nedge-1;j++) {
          //if ( !(edges[j].p1 < 0 && edges[j].p2 < 0) )
            for (int k=j+1;k<nedge;k++) {
              if ((edges[j].p1 == edges[k].p2) && (edges[j].p2 == edges[k].p1)) {
                edges[j].p1 = -1;
                edges[j].p2 = -1;
                
                edges[k].p1 = -1;
                edges[k].p2 = -1;
              }
              // Shouldn't need the following, see note above 
              if ((edges[j].p1 == edges[k].p1) && (edges[j].p2 == edges[k].p2)) {
                edges[j].p1 = -1;
                edges[j].p2 = -1;
                
                edges[k].p1 = -1;
                edges[k].p2 = -1;
              }
            }
        }
*/          
        
        /*
          Form new triangles for the current point
          Skipping over any tagged edges.
          All edges are arranged in clockwise order.
        */
        for (int j=0;j<nedge;j++) {
          if (edges[j].p1 == -1 || edges[j].p2 == -1)
            continue;
          if (ntri >= trimax) return -1;
          
          v[ntri].p1 = edges[j].p1;
          v[ntri].p2 = edges[j].p2;
          v[ntri].p3 = i;
          
          complete[ntri] = false;
          ntri++;
        }
      }
      // Include each point one at a time into the existing mesh (END)
      
      
      
      /*
        Final Step : Remove triangles with supertriangle vertices
        These are triangles which have a vertex number greater than nv
      */
      for (int i=0;i<ntri;i++) {
        if (v[i].p1 >= nv || v[i].p2 >= nv || v[i].p3 >= nv) {
          v[i] = v[ntri-1];
          ntri--;
          i--;
        }
      }
      INSTEAD - Just map filter them out...
      
      return ntri;
    }

/*    
    public static void main (String[] args) {
      int nv = 20;
      if (args.length > 0 && args[0] != null) nv = new Integer(args[0]).intValue();
      if (nv <= 0 || nv > 1000) nv = 20;
      
      //System.out.println("Creating " + nv + " random points.");
      
      XYZ[] points = new XYZ[ nv+3 ];
      
      for (int i=0; i<points.length; i++)
        points[i] = new XYZ( i*4.0, 400.0 * Math.random(), 0.0 );
      
      ITRIANGLE[]	 triangles 	= new ITRIANGLE[ nv*3 ];
      
      for (int i=0; i<triangles.length; i++)
        triangles[i] = new ITRIANGLE();
      
      int ntri = Triangulate( nv, points, triangles );
      
      // copy-paste the following output into free processing:
      //   http://processing.org/
      
      System.out.println("size(400,400); noFill();");
      
      for (int tt=0; tt<points.length; tt++){
        System.out.println("rect("+points[tt].x+","+points[tt].y+", 3, 3);");
      }
      
      for (int tt=0; tt<ntri; tt++){
        System.out.println("beginShape(TRIANGLES);");
        System.out.println("vertex( "+points[triangles[tt].p1].x+","+points[triangles[tt].p1].y+");");
        System.out.println("vertex( "+points[triangles[tt].p2].x+","+points[triangles[tt].p2].y+");");
        System.out.println("vertex( "+points[triangles[tt].p3].x+","+points[triangles[tt].p3].y+");");
        System.out.println("endShape();");
      }
      
    }
*/

  }
  
}

