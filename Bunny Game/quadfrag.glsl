#version 330 core

in vec4 fragWorldPos; // Received from the vertex shader

out vec4 FragColor;

uniform float scale; // You need to set this uniform to the desired scale for the checkerboard
uniform float offset; // You can set this uniform to offset the pattern if needed

void main() {
    // Calculate the checkerboard pattern
    bool x = int((fragWorldPos.x + offset) * scale) % 2 != 0;
    bool y = int((fragWorldPos.y + offset) * scale) % 2 != 0;
    bool z = int((fragWorldPos.z + offset) * scale) % 2 != 0;
    bool xorXY = x != y;
    bool checker = xorXY != z;

    // Set the color based on the checker pattern
    if (checker)
        FragColor = vec4(0, 0, 0, 1); // Black color
    else
        FragColor = vec4(1, 1, 1, 1); // White color
}
