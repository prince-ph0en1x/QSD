qubits 2

.init

h q0
h q1

#display

.quamdq(3)

ry q0,-0.927295
rz q0,-6.2832
rz q1,-3.141593

#displays

ry q1,0.927295
cnot q0,q1
ry q1,0.927295
cnot q0,q1

#display

rz q0,3.141593
rz q0,1.570796
rz q1,-1.570796
cnot q0,q1
rz q1,1.570796
cnot q0,q1
ry q0,-0.927295
rz q0,-3.141593

#display

h q0
h q1
x q0
x q1

h q0
cnot q1,q0
h q0

x q0
x q1
h q0
h q1

.showstate

display