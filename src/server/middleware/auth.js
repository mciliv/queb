const userMiddleware = (userService) => async (req, res, next) => {
  if (!userService) {
    return res.status(503).json({ error: "User service not available" });
  }

  const device_token = req.body.device_token || req.query.device_token;

  if (!device_token) {
    return res.status(400).json({ error: "Device token required" });
  }

  try {
    const user = await userService.getUserByDeviceToken(device_token);
    if (!user) {
      return res.status(404).json({ error: "User not found" });
    }
    req.user = user;
    next();
  } catch (error) {
    console.error("Authentication error:", error);
    res.status(500).json({ error: "Failed to authenticate user" });
  }
};

module.exports = { userMiddleware };



