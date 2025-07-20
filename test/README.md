
### **Run Simulations**

Use the `aldoclemente/fdapde-docker:latest` Docker image to compile and run C++ scripts.

#### 1. **Pull the Docker image**

```
docker pull aldoclemente/fdapde-docker:latest
```

#### 2. **Start the container**

```bash
./run_container
```

#### 3. **Navigate to the test directory inside the container**

```
cd ~/fdaPDE-testing/test
```

#### 4. **Compile and run a C++ script**

To compile and run a script:

```
./run_simulation.sh script.cpp
```

If your script requires **LBFGS++**, use:

```
./run_simulation_lbfgspp.sh script.cpp
```

