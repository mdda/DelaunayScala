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
  
  
}

