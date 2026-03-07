use ai_lib::prelude::*;

#[tokio::main]
async fn main() -> Result<(), AiLibError> {
    // Select your AI provider
    let client = AiClient::new(Provider::Groq)?;

    // Create a chat request
    let req = ChatCompletionRequest::new(
        client.default_chat_model(),
        vec![Message {
            role: Role::User,
            content: Content::Text("Explain transformers in one sentence.".to_string()),
            function_call: None,
        }]
    );

    // Send the request
    let resp = client.chat_completion(req).await?;

    // Get the response text
    println!("Answer: {}", resp.choices[0].message.content.as_text());
    Ok(())
}
