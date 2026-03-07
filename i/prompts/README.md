# LiteLLM Proxy Server

A unified gateway to 100+ LLM providers with the latest models pre-configured. Use this proxy from any language or application using the OpenAI SDK format.

## 🚀 Quick Start

### 1. Configure API Keys

Copy the example env file and add your API keys:

```bash
cp .env.example .env
# Edit .env with your actual API keys
```

### 2. Start the Proxy Server

```bash
docker compose up -d
```

The proxy will be available at: `http://localhost:4000`

### 3. Test It

```bash
curl http://localhost:4000/health
```

## 📋 Available Models

### Latest Models (Jan 2026)

| Model Name | Provider | Best For |
|------------|----------|----------|
| `gpt-4o` | OpenAI | General purpose, reasoning |
| `o1` | OpenAI | Complex reasoning tasks |
| `claude-3-7-sonnet` | Anthropic | Latest Claude, best performance |
| `claude-3-5-sonnet` | Anthropic | Balanced performance |
| `claude-3-5-haiku` | Anthropic | Fast, cost-effective |
| `gemini-2-0-flash` | Google | Fast multimodal |
| `gemini-1-5-pro` | Google | Long context (2M tokens) |
| `grok-2-vision` | xAI | Vision + reasoning |
| `deepseek-reasoner` | DeepSeek | Cost-effective reasoning |
| `deepseek-chat` | DeepSeek | Fast chat |

### Usage Aliases

- **Default chat**: `gpt-4o`
- **Reasoning**: `o1`
- **Vision**: `grok-2-vision`
- **Fast/cheap**: `gpt-4o-mini`

## 💻 Usage Examples

### Python

```python
from openai import OpenAI

client = OpenAI(
    base_url="http://localhost:4000",
    api_key="sk-1234"  # Your LITELLM_MASTER_KEY
)

response = client.chat.completions.create(
    model="claude-3-7-sonnet",
    messages=[{"role": "user", "content": "Hello!"}]
)

print(response.choices[0].message.content)
```

### Ruby

```ruby
require 'net/http'
require 'json'

uri = URI('http://localhost:4000/chat/completions')
request = Net::HTTP::Post.new(uri)
request['Authorization'] = 'Bearer sk-1234'
request['Content-Type'] = 'application/json'
request.body = {
  model: 'gpt-4o',
  messages: [{ role: 'user', content: 'Hello from Ruby!' }]
}.to_json

response = Net::HTTP.start(uri.hostname, uri.port) do |http|
  http.request(request)
end

puts JSON.parse(response.body)
```

### Node.js

```javascript
import OpenAI from 'openai';

const client = new OpenAI({
  baseURL: 'http://localhost:4000',
  apiKey: 'sk-1234'
});

const response = await client.chat.completions.create({
  model: 'claude-3-7-sonnet',
  messages: [{ role: 'user', content: 'Hello!' }]
});

console.log(response.choices[0].message.content);
```

### cURL

```bash
curl http://localhost:4000/chat/completions \
  -H "Authorization: Bearer sk-1234" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "gpt-4o",
    "messages": [{"role": "user", "content": "Hello!"}]
  }'
```

## 🔧 Configuration

### Add More Models

Edit `litellm_config.yaml`:

```yaml
model_list:
  - model_name: my-custom-model
    litellm_params:
      model: provider/model-name
      api_key: os.environ/MY_API_KEY
```

### Load Balancing

For the same model across multiple regions/providers:

```yaml
model_list:
  - model_name: gpt-4o
    litellm_params:
      model: azure/gpt-4o-eu
      api_base: https://eu.openai.azure.com/
  
  - model_name: gpt-4o
    litellm_params:
      model: azure/gpt-4o-us
      api_base: https://us.openai.azure.com/
```

### Environment Variables

Required:
- `LITELLM_MASTER_KEY` - Authentication key for your apps
- At least one provider API key (e.g., `OPENAI_API_KEY`)

Optional:
- `SLACK_WEBHOOK_URL` - Slack alerts for errors/budgets
- `AWS_*` - AWS Bedrock credentials
- `LITELLM_NUM_WORKERS` - Number of worker processes (default: 4)

## 📊 Monitoring

### Health Check

```bash
curl http://localhost:4000/health
```

### Prometheus Metrics

Available at: `http://localhost:4000/metrics`

### Grafana Dashboard

Open: `http://localhost:3000` (default user/pass: admin/admin)

### Logs

```bash
docker compose logs -f litellm
```

## 🛠️ Management

### Start/Stop

```bash
docker compose up -d      # Start
docker compose down       # Stop
docker compose restart    # Restart
```

### Update to Latest

```bash
docker compose pull
docker compose up -d
```

### View Models List

```bash
curl http://localhost:4000/models \
  -H "Authorization: Bearer sk-1234"
```

### Check Costs

```bash
curl http://localhost:4000/spend/tags \
  -H "Authorization: Bearer sk-1234"
```

## 🔒 Security

- **Change default keys**: Update `LITELLM_MASTER_KEY` in `.env`
- **Network isolation**: Proxy runs on localhost by default
- **API key encryption**: Keys encrypted at rest with `LITELLM_SALT_KEY`
- **Rate limiting**: Configure in `litellm_config.yaml`

## 🌐 Production Deployment

### Resource Requirements

- **Minimum**: 2 CPU, 4GB RAM
- **Recommended**: 4 CPU, 8GB RAM
- **High traffic**: 8+ CPU, 16GB RAM

### Deploy to Cloud

See [LiteLLM deployment docs](https://docs.litellm.ai/docs/proxy/deploy) for:
- Kubernetes/Helm
- AWS ECS
- Google Cloud Run
- Azure Container Apps

## 📚 API Documentation

### OpenAI-Compatible Endpoints

- `POST /chat/completions` - Chat completions
- `POST /completions` - Text completions
- `POST /embeddings` - Generate embeddings
- `POST /images/generations` - Generate images
- `GET /models` - List available models
- `GET /health` - Health check

Full API docs: `http://localhost:4000/docs`

## 🤝 Contributing

To add support for new models:

1. Check [LiteLLM providers](https://docs.litellm.ai/docs/providers)
2. Add to `litellm_config.yaml`
3. Add API key to `.env.example`
4. Update this README

## 📄 License

MIT

## 🔗 Links

- [LiteLLM Documentation](https://docs.litellm.ai/)
- [Supported Models](https://docs.litellm.ai/docs/providers)
- [GitHub Issues](https://github.com/BerriAI/litellm/issues)
