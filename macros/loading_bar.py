import sys

def loading_bar(event, tot_event, bar_length=50, process = "Vertexing"):
    progress = event/tot_event+1/tot_event
    arrow = '=' * int(bar_length * progress)
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write(f'\r[{arrow}{spaces}] {int(progress * 100)}%  {process}: event {event} of {int(tot_event)}   ')
    sys.stdout.flush()


def count_clusters(cl_list):
    n_cl = 0
    for cl in cl_list:
        if cl > 0:
            n_cl += 1
    return n_cl
