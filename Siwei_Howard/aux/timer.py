import time


class Timer:
    # We create a list to store every time period Timer recorded.
    results = []
    # We initialize the attributes t1 and t2 to be 0.
    # t1 records the time when Timer starts. t2 records the time when Timer ends.
    t1 = 0
    t2 = 0
    # We initialize the attribute duration to be 0.
    duration = 0
    # We default the unit to be 1. It could also be 1/60 if the unit is minutes, or 1/3600 if it's hours.
    unit = 1
    # We make the default unit_name to be 'seconds'.
    unit_name = 'seconds'

    # The constructor method:
    def __init__(self):
        pass

    # The start method.
    def start(self):
        # If t1 is 0, we start the Timer and record the current time in t1.
        if self.t1 == 0:
            self.t1 = time.time()
        # If t1 is not 0, the Timer is already started. We return error message.
        else:
            print('Error: the Timer has already started!')

    # The end method.
    def end(self):
        # If t1 is not 0, the Timer has already started. We end the Timer, and store the current time in t2.
        if self.t1 != 0:
            self.t2 = time.time()
            # We store the time duration between t1 and t2 in duration.
            self.duration = (self.t2 - self.t1) * self.unit
            # We print out the message of the time duration.
            #print('The time taken is ' + str(self.duration) + ' ' + self.unit_name + '.')
            # We append the time duration to the list.
            self.results.append(str(self.duration) + ' ' + self.unit_name)
            # After the Timer ends, we set t1, t2, and duration back to zero.
            self.t1 = 0
            self.t2 = 0
            self.duration = 0
        # If t1 is 0, the Timer is not currently running. We return an error message.
        else:
            print('Error: the Timer is not currently running!')

    # The method to configure the Timer to display either seconds, minutes, or hours.
    def set_timer(self, input):
        # If we take 's' as parameter, we set the unit of the Timer to be seconds.
        if input == 's':
            self.unit = 1
            self.unit_name = 'seconds'
        # If we take 'm' as parameter, we set the unit of the Timer to be minutes.
        if input == 'm':
            self.unit == 1/60
            self.unit_name = 'minutes'
        # If we take 'h' as parameter, we set the unit of the Timer to be hours.
        if input == 'h':
            self.unit = 1/3600
            self.unit_name = 'hours'

    # The method to retrieve the last timer result.
    def last_result(self):
        # We simply print out the last element of the list.
        return self.results[len(self.results) - 1]
