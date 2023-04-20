'''
@author: Karl J. Smith
'''
#TODO: Find a style guide for file headers. Then use style guide for file headers.

# Import libraries.
import time

#===============================================================================
# DEFINE SCRIPT FUNCTIONS
#===============================================================================

# Define function to update clock_dict with the function name and the time it takes to run the function.
def update_clock_dict(in_dict, func_name, func_start_time):
    
    if func_name not in in_dict.keys():
        in_dict[func_name] = [(time.time() - func_start_time)]
    else:
        in_dict[func_name].append((time.time() - func_start_time))

# Define function to print the clock_dict to the console.   
def print_clock_dict(in_dict):
    
    measured_time = sum([sum(in_dict[i]) for i in in_dict.keys()])
    
    def format_secs(secs):
        if type(secs) == int or type(secs) == float:
            if round(secs,2) > 0:
                return str(round(secs,2))
            else:
                return "< 0.01"
        else:
            return str(secs)
    
    fields_list = ["NAME","TIMES_CALLED","MAX_TIME","MIN_TIME","TOTAL_TIME","AVG_TIME","%_MEASURED"]
    print("{:<36} | {:<12} | {:<12} | {:<12} | {:<12} | {:<12} | {:<12}".format(*fields_list))
    
    for func in sorted(in_dict.keys()):
        
        times_called = len(in_dict[func])
        max_time = max(in_dict[func])
        min_time = min(in_dict[func])
        ttl_time = sum(in_dict[func])
        avg_time = ttl_time/times_called
        prcnt_total = (ttl_time/measured_time)*100
        
        row_list = [format_secs(i) for i in [func, times_called, max_time, min_time, ttl_time, avg_time, prcnt_total]]
        
        print("{:<36} | {:<12} | {:<12} | {:<12} | {:<12} | {:<12} | {:<12}".format(*row_list))

# Define function to print a 'percent complete' message to the console.
def print_prcnt_complete(pr_count, pr_total, increment):
    
    inc_list = [i/pr_total for i in range(0,100,increment)]
    
    if pr_count/pr_total in inc_list:
        
        print("Completed {}% of total.".format(str(int(pr_count/pr_total))))
        
    elif pr_count == pr_total:
        
        print("Completed 100% of total.")