import time
import thread  
class for_multi_thread(dict):
	def __ini__(self):
		print "HEHE"
	def job(self,input):
		for i in range(5):
			print "Threads",input,"time",i
			time.sleep(1)
		thread.exit_thread()
