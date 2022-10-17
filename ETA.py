
import time
import datetime
	
def ETA(on,of,t0,bigO):
	t1 = time.time()
	dt = t1-t0
	if (bigO == "N"):
		#linear time
		complete = (float(on)/float(of))
		time_per = dt/complete
		remaining = 1-complete
		time_remaining = time_per*remaining
		return (" working on "+str(on)+" of "+str(of)+", ETA based on O(N): "+str(datetime.timedelta(seconds=time_remaining)))
		
	if (bigO == "1/2N^2"):
		complete = (float(on)*(float(of)-float(on))+0.5*float(on)*float(on))
		time_per = dt/complete
		remaining = 0.5*float(of)*float(of)-complete
		time_remaining = time_per*remaining
		return (" working on "+str(on)+" of "+str(of)+", ETA based on O(1/2N^2): "+str(datetime.timedelta(seconds=time_remaining)))
		