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
    Triangulation subroutine
    Takes as input NV vertices in array pxyz
    Returned is a list of ntri triangular faces in the array v
    These triangles are arranged in a consistent clockwise order.
    The triangle array 'v' should be malloced to 3 * nv
    The vertex array pxyz must be big enough to hold 3 more points
    The vertex array must be sorted in increasing x values say
    
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

  //int Triangulate ( int nv, XYZ pxyz[], ITRIANGLE v[] ) {
  def Triangulation_bourke(measurements : List[Vector2]) : Seq[(Int, Int, Int)] = {
    case class ITRIANGLE(p1:Int, p2:Int, p3:Int)
    case class IEDGE(p1:Int = -1, p2:Int = -1)
    case class XYZ(p1:Float, p2:Float, p3:Float)

    val EPSILON = 0.000001

    /*
      Return TRUE if a point q(x,y) is inside the circumcircle made up of the points p1(x,y), p2(x,y), p3(x,y)
      The circumcircle centre (x,y) is returned and the radius r as a triple (x,y,r)
      NOTE: A point on the edge is inside the circumcircle
    */
    def CircumCircle( q:Vector2, p1:Vector2, p2:Vector2, p3:Vector2) : (Boolean, XYZ) {
      if ( Math.abs(p1.y-p2.y) < EPSILON && Math.abs(p2.y-p3.y) < EPSILON ) {
        System.out.println("CircumCircle: Points are colinear");
        (false, new XYZ(0,0,0))
      }
      else {
        val (mx1, my1) = ( (p1.x+p2.x)/2.0, (p1.y+p2.y)/2.0 )
        val (mx2, my2) = ( (p2.x+p3.x)/2.0, (p2.y+p3.y)/2.0 )
        
        val (m1,m2, xc,yc) = 
          if ( Math.abs(p2.y-p1.y) < EPSILON ) {
            val m2 = -(p3.x-p2.x) / (p3.y-p2.y)
            val xc =  (p2.x+p1.x) / 2.0
            val yc =  m2 * (xc - mx2) + my2
            (0,m2, xc,yc)
          }
          else 
            if ( Math.abs(p3.y-p2.y) < EPSILON ) {
              val m1 = -(p2.x-p1.x) / (p2.y-p1.y)
              val xc =  (p3.x + p2.x) / 2.0
              val yc =  m1 * (xc - mx1) + my1
              (m1,0, xc,yc)
            }
            else {
              val m1 = -(p2.x-p1.x) / (p2.y-p1.y)
              val m2 = -(p3.x-p2.x) / (p3.y-p2.y)
              val xc =  (m1*mx1 - m2*mx2 + my2 - my1) / (m1 - m2)
              val yc =  m1 * (xc - mx1) + my1
              (m1,m2, xc,yc)
            }
          
          val rsqr = {
            val (dx, dy) = (p2.x - xc, p2.y - yc)
            dx*dx + dy*dy
          }
          val drsqr = {
            val (dx, dy) = (q.x - xc, q.y - yc)
            dx*dx + dy*dy
          }
          
          ( drsqr <= rsqr, new XYZ(xc,yc, Math.sqrt(rsqr)) )
        }
      }

      // Got to here...
/*
      boolean complete[] 		= null;
      IEDGE 	edges[] 		= null;
      int 	nedge 			= 0;
      int 	trimax, emax 	= 200;
      int 	status 			= 0;
      
      boolean	inside;
      double 	xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r;
      double 	xmin, xmax, ymin, ymax, xmid, ymid;
      double 	dx, dy, dmax;
      
      int		ntri			= 0;
*/
      val nv = measurements.length
      val nedge = 0
      
      // Allocate memory for the completeness list, flag for each triangle
      trimax = 4*nv;
      complete = new boolean[trimax];
      for (int ic=0; ic<trimax; ic++) complete[ic] = false;
      
      /* Allocate memory for the edge list */
      edges = new IEDGE[emax];
      for (int ie=0; ie<emax; ie++) edges[ie] = new IEDGE();
      
      // Find the maximum and minimum vertex bounds.
      // This is to allow calculation of the bounding triangle
      val Pmin = Vector2( measurements.map(_.x).min, measurements.map(_.y).min )  // Top Left
      val Pmax = Vector2( measurements.map(_.x).max, measurements.map(_.y).max )  // Bottom Right
      val diameter = {
        val pdiff = Vector2( Pmax.x - Pmin.x, Pmax.y - Pmin.y)
        (pdiff.x > pdiff.y) ? pdiff.x : pdiff.y
      }
      val Pmid = Vector2( (Pmin.x + Pmax.x)/2.0, (Pmin.y + Pmax.y)/2.0 )
    
      /*
        Set up the supertriangle
        This is a triangle which encompasses all the sample points.
        The supertriangle coordinates are added to the end of the vertex list. 
        The supertriangle is the first triangle in the triangle list.
      */
      
      pxyz[nv+0].x = xmid - 2.0 * dmax;
      pxyz[nv+0].y = ymid - dmax;
      pxyz[nv+0].z = 0.0;
      
      pxyz[nv+1].x = xmid;
      pxyz[nv+1].y = ymid + 2.0 * dmax;
      pxyz[nv+1].z = 0.0;
      
      pxyz[nv+2].x = xmid + 2.0 * dmax;
      pxyz[nv+2].y = ymid - dmax;
      pxyz[nv+2].z = 0.0;
      
      v[0].p1 = nv;
      v[0].p2 = nv+1;
      v[0].p3 = nv+2;
      
      complete[0] = false;
      ntri = 1;
      
      // Include each point one at a time into the existing mesh
      for (int i=0;i<nv;i++) {
        xp = pxyz[i].x;
        yp = pxyz[i].y;
        nedge = 0;
        
        /*
          Set up the edge buffer.
          If the point (xp,yp) lies inside the circumcircle then the
          three edges of that triangle are added to the edge buffer
          and that triangle is removed.
        */
        XYZ circle = new XYZ();
        for (int j=0;j<ntri;j++) {
          if (complete[j])
            continue;
          x1 = pxyz[v[j].p1].x;
          y1 = pxyz[v[j].p1].y;
          
          x2 = pxyz[v[j].p2].x;
          y2 = pxyz[v[j].p2].y;
          
          x3 = pxyz[v[j].p3].x;
          y3 = pxyz[v[j].p3].y;
          
          inside = CircumCircle( xp, yp,  x1, y1,  x2, y2,  x3, y3,  circle );
          
          xc = circle.x; yc = circle.y; r = circle.z;
          
          if (xc + r < xp) 
            complete[j] = true;
            
          if (inside) {
            /* Check that we haven't exceeded the edge list size */
            if (nedge+3 >= emax) {
              emax += 100;
              IEDGE[] edges_n = new IEDGE[emax];
              for (int ie=0; ie<emax; ie++) edges_n[ie] = new IEDGE();
              System.arraycopy(edges, 0, edges_n, 0, edges.length);
              edges = edges_n;
            }
            edges[nedge+0].p1 = v[j].p1;
            edges[nedge+0].p2 = v[j].p2;
            
            edges[nedge+1].p1 = v[j].p2;
            edges[nedge+1].p2 = v[j].p3;
            
            edges[nedge+2].p1 = v[j].p3;
            edges[nedge+2].p2 = v[j].p1;
            
            nedge += 3;
            
            v[j].p1 = v[ntri-1].p1;
            v[j].p2 = v[ntri-1].p2;
            v[j].p3 = v[ntri-1].p3;
            
            complete[j] = complete[ntri-1];
            
            ntri--;
            j--;
          }
        }
        
        /*
          Tag multiple edges
          Note: if all triangles are specified anticlockwise then all
          interior edges are opposite pointing in direction.
        */
        for (int j=0;j<nedge-1;j++) {
          //if ( !(edges[j].p1 < 0 && edges[j].p2 < 0) )
            for (int k=j+1;k<nedge;k++) {
              if ((edges[j].p1 == edges[k].p2) && (edges[j].p2 == edges[k].p1)) {
                edges[j].p1 = -1;
                edges[j].p2 = -1;
                edges[k].p1 = -1;
                edges[k].p2 = -1;
              }
              /* Shouldn't need the following, see note above */
              if ((edges[j].p1 == edges[k].p1) && (edges[j].p2 == edges[k].p2)) {
                edges[j].p1 = -1;
                edges[j].p2 = -1;
                edges[k].p1 = -1;
                edges[k].p2 = -1;
              }
            }
        }
        
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
      
      
      /*
        Remove triangles with supertriangle vertices
        These are triangles which have a vertex number greater than nv
      */
      for (int i=0;i<ntri;i++) {
        if (v[i].p1 >= nv || v[i].p2 >= nv || v[i].p3 >= nv) {
          v[i] = v[ntri-1];
          ntri--;
          i--;
        }
      }
      
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

