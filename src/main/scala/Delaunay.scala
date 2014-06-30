//package Delaunay

import VSFM.{Vector2,Vector3}

package object Delaunay {
  def Triangulation(measurements : List[Vector2]) : Seq[(Int, Int, Int)] = {
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
}


/*
//import java.awt.Point;
//import java.util.Iterator;
//import java.util.TreeSet;

// http://stackoverflow.com/questions/5825089/how-does-this-code-for-delaunay-triangulation-work
class DelaunayTriangulation(size: Int) {
   int[][] adjMatrix;

   DelaunayTriangulation(int size) {
     this.adjMatrix = new int[size][size];
   }
   public int[][] getAdj() {
     return this.adjMatrix;
   }

   def getEdges(int n, int[] x, int[] y, int[] z) : Set {
     if (n == 2) {
       this.adjMatrix[0][1] = 1;
       this.adjMatrix[1][0] = 1;
       result.add(new GraphEdge(new GraphPoint(x[0], y[0]), new GraphPoint(x[1], y[1])));

       return result;
     }
     
     TreeSet result = new TreeSet();

     for (int i = 0; i < n - 2; i++) {
       for (int j = i + 1; j < n; j++) {
         for (int k = i + 1; k < n; k++) {
           if (j == k) {
             continue;
           }
           
           int xn = (y[j] - y[i]) * (z[k] - z[i]) - (y[k] - y[i]) * (z[j] - z[i]);
           int yn = (x[k] - x[i]) * (z[j] - z[i]) - (x[j] - x[i]) * (z[k] - z[i]);
           int zn = (x[j] - x[i]) * (y[k] - y[i]) - (x[k] - x[i]) * (y[j] - y[i]);
           
           boolean flag;
           if (flag = (zn < 0 ? 1 : 0) != 0) {
             for (int m = 0; m < n; m++) {
               flag = (flag) && ((x[m] - x[i]) * xn + (y[m] - y[i]) * yn + (z[m] - z[i]) * zn <= 0);
             }

           }

           if (!flag) {
             continue;
           }
           
           //System.out.println(x[i]+" "+ y[i] +"----"+x[j]+" "+y[j]);
           result.add(new GraphEdge(new GraphPoint(x[i], y[i]), new GraphPoint(x[j], y[j])));
           this.adjMatrix[i][j] = 1;
           this.adjMatrix[j][i] = 1;
           
           //System.out.println(x[j]+" "+ y[j] +"----"+x[k]+" "+y[k]);
           result.add(new GraphEdge(new GraphPoint(x[j], y[j]), new GraphPoint(x[k], y[k])));
           this.adjMatrix[j][k] = 1;
           this.adjMatrix[k][j] = 1;
           
           //System.out.println(x[k]+" "+ y[k] +"----"+x[i]+" "+y[i]);
           result.add(new GraphEdge(new GraphPoint(x[k], y[k]), new GraphPoint(x[i], y[i])));
           this.adjMatrix[k][i] = 1;
           this.adjMatrix[i][k] = 1;
           
           //System.out.println("----------");
         }

       }

     }

     return result;
   }
*/
/*
   public TreeSet getEdges(TreeSet pointsSet) {
     if ((pointsSet != null) && (pointsSet.size() > 0)) {
       int n = pointsSet.size();

       int[] x = new int[n];
       int[] y = new int[n];
       int[] z = new int[n];

       int i = 0;

       Iterator iterator = pointsSet.iterator();
       while (iterator.hasNext()) {
         Point point = (Point)iterator.next();

         x[i] = (int)point.getX();
         y[i] = (int)point.getY();
         z[i] = (x[i] * x[i] + y[i] * y[i]);

         i++;
       }

       return getEdges(n, x, y, z);
     }

     return null;
   }
}
*/   
