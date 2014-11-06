
%simple demo script for matlab dce
%test by Mike Griffiths 22nd January 2007
resource = findResource('scheduler', 'configuration', 'sge');
set(resource, 'configuration', 'sge');

parjob=createParallelJob(resource);

set(parjob,'MinimumNumberOfWorkers',numprocs);
set(parjob,'MaximumNumberOfWorkers',numprocs);
%createTask(parjob, 'rand', 1, {3});
%createTask(parjob, 'myrandpar', 2, {5});
createTask(parjob, 'mygetpictfunc');

tic
submit(parjob);

waitForState(parjob);
myruntime=toc

