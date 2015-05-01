'''
iPy contains helpers for iPython parallel computing tools
'''
import os,sys,time


def run_taskqueue(tasks,tc=None,sleeptime=60):
    '''
    run_taskqueue takes a list of IPython client.ITask objects (eg client.StringTask or client.MapTask),
    a TaskClient (or if not specified attempts to construct the default IPython.kernel.client.TaskController())
    and an optional sleeptime=seconds_to_sleep

    returns a list of completed task_result objects
    '''
    
    if not tc:
        from IPython.kernel import client
        tc = client.TaskController()

    tc.clear()

    taskids = []
    for t in tasks:
        taskids.append(tc.run(t))
    print taskids

    #left = len([ tc.get_task_result(taskid) for taskid in taskids  if tc.get_task_result(taskid) is None])
    left = tc.queue_status()['pending']+ tc.queue_status()['scheduled']

    while left > 0:
    #left = len([ tc.get_task_result(taskid) for taskid in taskids  if tc.get_task_result(taskid) is None])
        print >> sys.stderr, tc.queue_status()
        time.sleep(sleeptime)
        #left = len([ tc.get_task_result(taskid) for taskid in taskids  if tc.get_task_result(taskid) is None])
        left = tc.queue_status()['pending']+ tc.queue_status()['scheduled']

    return [ tc.get_task_result(taskid) for taskid in taskids ]
