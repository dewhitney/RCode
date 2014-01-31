log.entry = function(hrs,min,sec,miles)
{
  t = hrs*60*60 + min*60 + sec
  s = t/miles
  p = paste(round(s)%/%60,":",round(s)%%60, sep="")
  time = paste(hrs,min,sec, sep=":")
  entry = c(dist=miles,time=time,pace=p)
  return(entry)
}