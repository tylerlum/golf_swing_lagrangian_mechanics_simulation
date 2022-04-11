# golf_swing_lagrangian_mechanics_simulation

The purpose of this repository is to create a numerical simulation of a golf swing with a Lagrangian mechanics model.

## Visualizations of Golf Swings

### Passive Arm System: Golfer Simply Allows Gravity to Create the Swing

![passive](https://user-images.githubusercontent.com/26510814/162670542-5f607c12-73d3-4930-8198-e1d0e6e337b6.gif)

### Controlled Arm System: Low Arm Angular Acceleration

![controlled_a2](https://user-images.githubusercontent.com/26510814/162670539-3cadaea2-6faf-4569-9ba3-6cd518259d9e.gif)

### Controlled Arm System: Medium Arm Angular Acceleration

![controlled_a5](https://user-images.githubusercontent.com/26510814/162670530-2574daf9-16ba-4956-bd45-108efca53a2b.gif)

### Controlled Arm System: High Arm Angular Acceleration

![controlled_a20](https://user-images.githubusercontent.com/26510814/162670536-add72abf-af55-4ce4-ba36-4094b9357c03.gif)

## Optimization of Wrist Flick Angle

In this model, we initially have a fixed wrist angle of 90 degrees, and then allow the wrist angle to move freely at an arm angle of $\theta_{fixed}$. We optimize this angle so maximize the clubhead speed. We find that for higher arm angular accelerations, the optimal wrist flick angle is lower (later time).

### Passive Arm System: Golfer Simply Allows Gravity to Create the Swing

![clubhead_speed_vs_theta_fixed](https://user-images.githubusercontent.com/26510814/162670989-691c0725-ef2e-4d28-8fce-68ae22a8d794.png)

### Controlled Arm System: Low Arm Angular Acceleration

![clubhead_speed_vs_theta_fixed_a_2](https://user-images.githubusercontent.com/26510814/162670992-a08d396f-ab27-45d4-ac72-22497c9f3747.png)

### Controlled Arm System: Medium Arm Angular Acceleration

![clubhead_speed_vs_theta_fixed_a_5](https://user-images.githubusercontent.com/26510814/162670983-59be6a58-0e51-46ba-9f7a-cde44a6d78f6.png)

### Controlled Arm System: High Arm Angular Acceleration

![clubhead_speed_vs_theta_fixed_a_20](https://user-images.githubusercontent.com/26510814/162670990-19b96775-e466-46af-a059-a692bfcb216c.png)

## Progression of Arm and Wrist Angles

The progression of the arm and wrist angles for various systems appear to have the same general shape.

### Passive Arm System: Golfer Simply Allows Gravity to Create the Swing

![passive](https://user-images.githubusercontent.com/26510814/162670942-e36980b3-3df0-4d7d-8925-f8a29d092e02.png)

### Controlled Arm System: Low Arm Angular Acceleration

![a2](https://user-images.githubusercontent.com/26510814/162670943-22c4b37a-51f7-4bba-afc4-175409fe452a.png)

### Controlled Arm System: Medium Arm Angular Acceleration

![a5](https://user-images.githubusercontent.com/26510814/162670945-6a4a5825-2fe2-473f-add0-54050cb2f91a.png)

### Controlled Arm System: High Arm Angular Acceleration

![a20](https://user-images.githubusercontent.com/26510814/162670944-4ee389e8-b84d-4f42-a88a-13dbd19b412d.png)


