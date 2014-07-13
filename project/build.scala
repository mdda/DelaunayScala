import sbt._
import Keys._

object DelaunayScalaBuild extends Build {
  def scalaSettings = Seq(
    scalaVersion := "2.10.0",
    scalacOptions ++= Seq(
      "-optimize",
      "-unchecked",
      "-deprecation"
    )
    
  )
  
  def buildSettings =
    Project.defaultSettings ++
    scalaSettings

  lazy val root = {
    val settings = buildSettings ++ Seq(name := "DelaunayScala")
    Project(id = "DelaunayScala", base = file("."), settings = settings)
  }
  
}
