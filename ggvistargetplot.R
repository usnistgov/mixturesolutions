###make the plot with ggvis so it can be interactive

load("exampleM1.rData")
g%>%ggvis(~value, ~value2,fill=~factor(variable))%>%
  layer_points()%>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.03,npoints=25)),~x,~y,fill:="green",opacity:=0.5)%>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.09,npoints=25)),~x,~y,fill:="red",opacity:=0.3) %>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.06,npoints=25)),~x,~y,fill:="yellow",opacity:=0.3)
stp2%>%ggvis(~value, ~value2,fill=~factor(variable))%>%
  layer_points()%>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.03,npoints=25)),~x,~y,fill:="green",opacity:=0.5)%>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.09,npoints=25)),~x,~y,fill:="red",opacity:=0.3) %>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.06,npoints=25)),~x,~y,fill:="yellow",opacity:=0.3)


stp2%>%ggvis(~value, ~value2,fill=~factor(variable))%>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.03,npoints=25)),~x,~y,fill:="green",opacity:=0.5)%>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.09,npoints=25)),~x,~y,fill:="red",opacity:=0.3) %>%
  layer_paths(data = data.frame(circleFun(center=c(0,0),diameter=0.06,npoints=25)),~x,~y,fill:="yellow",opacity:=0.3)
  layer_points(key=~id)%>%
  add_tooltip(all_values,"hover")
  #the tooltip is still not working because add tooltip needs to refer to a vis and i dont know how ; look it up
all_values<-function(x){if(is.null(x))return(NULL)
            row<-stp2[stp2$id==x$id ,]
            paste0(names(row), ": ", format(row), collapse="<br />")}
#this function is supposed to provide the hover-over of all of the values, but currently
#isn't getting the ID variable returned properly.